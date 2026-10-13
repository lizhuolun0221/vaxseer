# VaxSeer 复现环境

环境名称：`vaxseer_lzl`。代码基线：`470a25b`。

## 激活

```bash
source "/home/lizhuolun/Nat Med/vaxseer_reproduction/activate_vaxseer_lzl.sh"
conda env list
python --version
```

环境位于项目内的 `vaxseer_lzl/`，是独立 Conda 环境。当前执行权限不能写入 `~/.conda/environments.txt`，所以默认终端的环境列表可能看不到它。激活脚本为当前终端设置 `CONDA_ENVS_PATH`，使 Conda 能按名称发现环境，同时修正此机器上系统 Python 优先于激活环境的 PATH 顺序。若在 IDE 中选择解释器，使用：

`/home/lizhuolun/Nat Med/vaxseer_reproduction/vaxseer_lzl/bin/python`

## 安装记录

- Python 3.8.18、pip 24.0，通过 conda-forge 安装。
- 官方 environment.yaml 的 Python 依赖版本均保留；科学计算库通过 pip 安装相同版本。
- PyTorch 1.13.1+cu117、Torchvision 0.14.1+cu117、Torchaudio 0.13.1+cu117：从上海交大镜像下载，与 PyTorch 官方索引提供的 SHA256 完全一致。
- 其他 pip 依赖使用清华镜像安装。
- 固定 NumPy 1.24.4；pip 24.0 兼容旧版 Lightning 的依赖元数据，不要直接升级 pip。
- 固定 Conda libsqlite 3.45.3，解决旧 PyTorch 与新 SQLite/ICU 同进程导入时的 C++ ABI 冲突。
- MMseqs2 官方 Linux AVX2 二进制，版本 `a13af93468c997d73572c9326b69ac0c66c6583a`，安装在环境 `bin/mmseqs`。需要支持 AVX2 的 CPU。
- 安装日志和下载缓存保存在 `.setup/`、`.conda-pkgs/`，与环境目录一起排除在 Git 状态之外。

`conda-explicit.txt` 记录 Conda 软件包；`requirements-lock.txt` 记录实际 Python 软件包；`requirements.txt` 记录顶层依赖。MMseqs2 是单独安装的二进制，不包含在这些包清单中。

## 验证

在仓库根目录执行：

```bash
source activate_vaxseer_lzl.sh
python environment_setup/verify_environment.py
cd vaxseer
python -m bin.train --help
```

验证结果保存在 `verification.json`：19 个库/项目模块导入、CPU 反向传播、pip check、MMseqs2 版本命令均通过。训练入口 `--help` 退出码为 0。

PyTorch 编译目标为 CUDA 11.7，但当前会话 `torch.cuda.is_available()` 为 false，`nvidia-smi` 无法连接驱动。尚未验证 GPU 训练，也未下载完整数据或模型权重，未执行论文结果复现。
