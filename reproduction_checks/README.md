# VaxSeer 运行检查

本目录保存检查脚本、日志、配置和输出。所有命令在 `vaxseer_reproduction/` 执行，先激活已有环境：

```bash
source ./activate_vaxseer_lzl.sh
```

## 严格加载与 CPU/GPU 前向

```bash
python reproduction_checks/check_models.py \
  --subtype a_h1n1 --year 2018 --gpu 1 \
  --output reproduction_checks/h1n1_2018_CHECK_NEW
```

已通过结果：`h1n1_2018_20260917/summary.json`。两类检查点按官方类严格加载，使用官方 H1N1 参考 HA（566 aa，含 2 个 X）完成 FP32 CPU/GPU 前向。dominance 检查时间 bin 90–92；antigenicity 使用已天然对齐的同一参考序列配对及官方 dataset/collator。最大 CPU/GPU 差异分别为 `2.2531e-5`、`1.1921e-7`，均通过所记录的容差。

运行脚本设为离线，仅使用已有权重和缓存配置。首次交互检查由 Transformers 补取过 GPT-2 配置（665 字节），没有补取模型权重。检查不训练、不生成序列、不修改检查点；它验证可执行性，不验证预测准确性。

## 小规模完整链路

```bash
python reproduction_checks/run_smoke.py \
  --dataset ../Full_data/datasets/h1n1_standard_20260917 \
  --output reproduction_checks/h1n1_2018_SMOKE_NEW \
  --year 2018 --gpu 1 --size 4
```

已通过结果：`h1n1_2018_smoke_20260917/summary.json`。按 metadata 的提交日期排除 2018-02-01 及以后、以及无完整提交日期的记录；然后使用官方预处理，采样时间窗为 `[2015-02, 2018-02)`，人源、长度至少 553、出现次数至少 5、bin 大小至少 1000。保留 metadata 部分日期，官方预处理排除只有年/月的采样日期。

预处理得到 214 条不同 HA 序列。本次仅按输出顺序取前 4 条，不根据模型预测筛选。使用 MMseqs2 对齐和官方模型推理入口，再使用官方评分脚本；生成 16 个配对、12 行 dominance 预测（4 条 × 3 个时间 bin）、16 行 antigenicity 预测及 4 行评分。输出均有限，每条评分的病毒概率支持总和为 1，并通过独立加权重算核对。

官方 `build_test_pairs.py` 默认排除相同 query/target 编号。本地新增可选 `--include_self_pairs`，默认行为保持不变。小规模检查启用该参数，以保留完整配对矩阵；另存 `pairs_official_default.csv` 验证原行为为 12 对，仅少 4 个对角配对。使用 `run_smoke.py --official-exclude-self` 可检查官方设置。正式比较必须记录并统一此设置，不能把两种设置的结果混合。

官方 score 是以预测流行权重加权的预测抗原距离（训练标签为 log2 HI 降低量），数值较小代表模型预测的抗原差异较小，不是百分比覆盖率。

本次是执行检查，不是完整历史回顾评估。标准数据没有限定为论文编号集；当前 metadata 的历史提交日期也不能恢复历史数据库版本。尚未核对真实 HI 配对覆盖、HI 报告时间、原模型测试划分、未来流行权重、WHO 候选和基线。不得将这 4 条的结果当作论文指标或候选推荐。

## 路径适配

新的检查脚本按真实目录定位唯一检查点，并检查配对矩阵。官方旧入口存在目录/文件名不一致：dominance 示例路径缺少 `weight_loss_by_count/`，预处理当前输出含 `minLen553`，与旧示例的 `lenQuantile0.2` 不同。新的脚本使用真实输出路径，避免这些不一致；不运行旧脚本里的日志删除操作。

输出目录必须使用新名称。`--resume-preprocessing` 仅用于已经准备好 cohort 和预处理文件、还未生成配对的目录，不能覆盖已完成的推理。

## HI 名称映射、标签与真实输入（2026-09-18）

已完成证据见 `hi_h1n1_2018_20260918/CHECK_REPORT.md`：29,203 个官方兼容标签逐对核对通过，6,904 个唯一映射的历史提交配对全部对齐，64 对真实 HI 的 GPU 输入检查通过。原训练/测试划分未恢复，没有报告泛化指标。

重建命令（输出须使用新名称）：

```bash
python reproduction_checks/prepare_hi_audit.py \
  --dataset ../Full_data/datasets/h1n1_standard_20260917 \
  --cutoff 2018-02-01 --paper-only \
  --output reproduction_checks/hi_h1n1_2018_NEW

python reproduction_checks/check_hi_inputs.py \
  --audit reproduction_checks/hi_h1n1_2018_NEW \
  --dataset ../Full_data/datasets/h1n1_standard_20260917 \
  --run-prediction --gpu 1

python reproduction_checks/triage_hi_missing.py \
  --audit reproduction_checks/hi_h1n1_2018_NEW \
  --sequence-index ../Full_data/catalog/audit_20260916/sequences.sqlite \
  --metadata-index ../Full_data/catalog/meta_h1n1_20260916/metadata.sqlite
```

仅采样日期条件保存官方编号展开加权标签；提交日期条件同时保存官方兼容标签和唯一映射、源 CSV 行等权的保守标签。两套标签不得混用。逐名称/逐源行/逐序列对的来源表、缺失清单和 notebook 均保留在输出目录。原始同源 HI 值缺失的两种参考名称影响 49 行，已隔离。

缺失分类清单中，228 个已有论文 HA 只需优先补 metadata；`missing_triage/` 分别保存 EPI 和 EPI_ISL 列表，不能混合使用。另有 1,828 个待检索名称，不代表 GISAID 一定没有这些记录。
