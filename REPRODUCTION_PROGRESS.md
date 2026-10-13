# VaxSeer 复现进度（更新：2026-09-18）

## 已完成

- H1N1-2018 dominance、antigenicity 实际严格加载、完整参考 HA 的 CPU/GPU 前向和结果一致性检查。
- 标准 H1N1 HA/metadata 导出：178,735 个编号一一对应，官方解析器验证通过；保留部分日期，隔离 3 个日期冲突。完整原始记录保持不变。
- 使用 4 条历史提交日期筛选后的 HA 跑通 MMseqs2 → dominance → antigenicity → 官方 coverage score。16 个配对、12 行进化预测、16 行抗原性预测、4 行评分；独立重算通过。
- 两个可重复运行的模型检查脚本、一个标准数据导出脚本及配置/来源/结果文件。
- 修复原始序列审计的发现范围，排除 `catalog/`、`datasets/`、`tools/`、`download_lists/`，防止下次把派生 FASTA 当作新增下载；验证仍只发现原始 24 个 FASTA。
- 配对脚本新增保留相同编号配对的可选开关；验证默认 12 对、新设置 16 对。官方默认行为保留。
- H1N1 HI 名称映射、缺失/多序列歧义清单和源行到蛋白序列对的来源记录已生成；本轮限定论文 HA 编号，分别检查仅采样日期及同时核对提交日期的两套条件。
- 29,203 个官方兼容标签已与实际运行的原版 `prepare.py` 按蛋白哈希核对，配对集合完全一致，标签最大差异约 `2e-14`。
- 2018-02 前同时满足采样/提交日期条件、两端唯一映射的 6,904 个保守 HI 配对已全部对齐（1,178 条 virus 蛋白、24 条 reference 蛋白，缺配 0）；64 个真实 HI 输入完成 GPU 推理检查，输出有限且标签不变。该保守协议对源 CSV 行等权，与官方编号多重性加权分开保存。
- 已定位优先补 metadata 的 228 个本地已有论文 HA 编号/228 个毒株编号（涉及 136 个 HI 名称）；另有 1,828 个 HI 名称未匹配到本地长度合格 HA，已保存检索清单。原始 FASTA 仍为已索引的 24 个，无新增或修改。
- HI 复核 notebook 已执行成功，公式重算及配对唯一性检查通过。

运行证据见 `vaxseer_reproduction/reproduction_checks/README.md` 及两个运行目录的 `summary.json`；数据检查见 `Full_data/datasets/h1n1_standard_20260917/CHECK_REPORT.md`。

本轮 HI 证据见 `vaxseer_reproduction/reproduction_checks/hi_h1n1_2018_20260918/CHECK_REPORT.md`；完整对齐输入为该目录下 `input_checks/aligned_unique_sequence_pairs.csv`；已执行 notebook 为 `hi_match_review.executed.ipynb`。


## 近期补充：dominance 回顾统计与论文频率指标对照（2026-09-20）

- 完成 H1N1-2018 执行检查的预测结果归一化、按 3 个两月预测 bin 分别排序，并以 2018-10 至 2019-03 标准数据作 HA 完整蛋白序列精确匹配。214 条历史候选在三段中的匹配记录数分别为 66、53、17；对应观察到的不同候选 HA 数为 8、7、5。候选集外样本较多，匹配记录占同期全部合格样本比例分别为 3.72%、0.95%、0.44%。统计脚本及结果在 `vaxseer_reproduction/reproduction_checks/h1n1_2018_all214/observed_comparison/`。
- `*_fixed214.csv` 固定同一批 214 条候选，以未来 HA 氨基酸序列完全一致为匹配规则，未设置未来出现次数阈值，零次候选保留。`observed_percent_within_214` 是在 214 条匹配记录内部的占比；它不等于全体未来样本中的流行频率。原 `*_comparison.csv` 的 `observed_percent_all_eligible` 使用同期全部合格记录作分母。
- 对照论文的 dominance 频率评估：仓库官方 `evaluation/jupyter_notebooks_figs/Appendix_eval_freq.ipynb` 将预测输入限制为截止前 6 个月的 HA，在候选集合中按多个预测 bin 的概率聚合排序；对真实未来季节序列按完整蛋白序列精确匹配，未匹配记 0；计算预测前 20% 序列在未来季节全部合格 HA 样本中的频率总和（top-20% frequency coverage）。该 notebook 当前代码示例仅配置 H3N2、2014–2019；论文正文报告 dominance 方法比较及单侧 Wilcoxon 检验。
- 因而当前 214 条工作与原论文指标是**同一类回顾性频率问题，但统计量与设置不相同**：当前使用三年候选集、分开两月窗口、报告单序列计数/排名；原 notebook 使用前 6 个月序列、聚合预测、top 20% 总频率且以全部未来样本为分母。当前结果不能直接据此复述论文的 top-20% 结果。
- 下一步如需同口径比较，应重建截止点前 6 个月的历史输入候选，对每个候选聚合预测季节 bin 的概率，按模型排序取前 20%，用目标季节（10 月至次年 3 月）的全体合格数据计算这些序列的真实频率和；与 `Last`（历史频率）等基线使用相同候选支持、未来分母及年份，跨年份报告结果。当前保守的 metadata 提交日期限制、数据覆盖不全和历史快照不一致仍需标明。

## 下一阶段

1. 核对 HI 来源报告时间和原模型评估划分，再构建正式评估输入。不能凭重建的随机划分声称恢复原测试集；当前官方配对 CSV 只有 virus/reference/hi，不能直接证明每个 HI 值在历史截止日前已发表。现有 6,904 个历史重建配对可能与模型训练数据重叠，暂不作为独立测试集。
2. 整理 WHO 候选、未来流行时窗和基线共享输入；固定官方标签或保守等权标签协议。2018 提交日期条件中，官方编号多重性加权影响 714 个序列对；不可在不同方法之间混用标签或配对支持。
3. 补齐 HI 关键 metadata：优先使用 `vaxseer_reproduction/reproduction_checks/hi_h1n1_2018_20260918/missing_triage/local_paper_isolates_needing_metadata.txt`（228 个 EPI_ISL），HA 编号另存 `local_paper_ha_needing_metadata.txt`（228 个 EPI），不能混合上传。此批已有 FASTA，不必重复下载序列。
4. 用相同候选集、未来流行序列、配对支持、时间划分和指标评估预训练模型、基线与 WHO 历史候选，保存样本数和不确定性；然后扩展年份。
5. 继续补论文序列和 metadata。当前原始文件 inventory 与已有索引一致，论文 HA 编号仍覆盖 210,070/394,075（53.3071%），缺 184,005；H3N2 metadata 尚未提供，H1N1 metadata 缺 FASTA 10,120 个、FASTA 缺 metadata 21,363 个。
6. H1N1-2012 dominance 权重继续按约定跳过。

本轮得到的是已验证的小规模执行闭环，不是论文指标复现。2018 小规模使用当前记录中采样和提交日期符合要求的部分数据；这不足以恢复当时数据库版本，也未将全部标准输入限定为论文编号集。

2026-09-18 的 HI 重建已限定论文编号，但仍缺报告日期和原测试划分。真实 HI 的 64 对检查是刻意平衡的输入检查，不能报告为模型泛化指标。

## 尚待适配的官方旧脚本

- `evaluation/pipeline/run_vaxseer.sh` dominance 路径缺 `weight_loss_by_count/`，并使用与当前预处理输出不同的文件名；本轮新脚本已按真实路径适配。
- `evaluation/pipeline/run_sweep_comb.sh` 使用 `seq 2024 2021`，默认步长下为空；批量正式评估前需更正年份范围，并处理明确跳过的模型。
- coverage 脚本会对实际存在 HI 配对的病毒概率支持重新归一化。正式比较必须核对不同候选/方法的支持集合，不能忽略缺配。
