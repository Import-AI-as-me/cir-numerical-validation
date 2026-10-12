# CIR 定价与敏感度数值验证

**价格完全匹配，敏感度仍可能出错。** 本项目的解析反例中，连续时间 CIR 债券价格与年度网格近似都等于 **0.8198175474**，但网格对初始短率的敏感度高估 **3.06237%**。

这是一项可复现的数值验证研究：分开检查积分偏差、状态离散偏差和 Monte Carlo 采样误差，并提供原始结果与实现检查。项目复现已有数值机制，不宣称新定理或交易盈利。

[English](README.md) · [英文技术报告](report.md) · [中文研究报告 PDF](report.zh-CN.pdf)

![年度网格在构造的初始短率处匹配价格，但初始短率敏感度仍高估约 3.06%。](price_greek_cancellation.png)

## 核心例子

考虑 CIR 短率模型，价格为折现因子的期望，敏感度定义为初始短率导数的相反数：

$$
dr_t=\kappa(\theta-r_t)\,dt+\sigma\sqrt{r_t}\,dW_t,
\qquad P=\mathbb E\!\left[e^{-\int_0^T r_s\,ds}\right],
\qquad G=-\frac{\partial P}{\partial r_0}.
$$

取 $\kappa=0.6$、$\theta=0.04$、$\sigma=0.1$、期限 $T=5$ 年，使用年度网格。在人为构造的初值 $r_0^*=0.0400853469$ 处：

| 输出 | 连续时间真值 | 精确网格转移＋梯形积分 |
| --- | ---: | ---: |
| 单位本金价格 $P$ | 0.8198175474 | 0.8198175474 |
| 初始短率敏感度 $G$ | 1.2852865306 | 1.3246467181 |
| 敏感度相对误差 | — | +3.06237% |

价格相等来自仿射系数的抵消；保存的双精度计算得到零价格绝对误差。这是验证用反例，不代表其在市场中的发生频率。精确模拟离散时点的状态，也没有消除路径积分误差。推导见[技术报告](report.md)，完整数值见[反例表](cancellation_counterexamples.csv)。

## 本地运行

使用 **Python 3.12、NumPy 和 Matplotlib**。数据全部为合成数据；安装依赖后不需要行情订阅、API 密钥或网络连接。

```bash
git clone https://github.com/Import-AI-as-me/cir-numerical-validation.git
cd cir-numerical-validation
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

运行实现测试与数值检查（`--verify-only` 只将检查结果以 JSON 打印到终端，不写结果文件）：

```bash
python -m unittest -v test_validation
python study.py --verify-only
```

[验证记录](validation.json)：在新的 Python 3.12 环境中，9 项实现测试、数值检查及包含 36 个组合的快速运行均已通过。

小规模检查运行流程：

```bash
python study.py --quick --output-dir results/quick
```

复现完整实验配置：

```bash
python study.py --reps 80 --paths 8000 --output-dir results/full
```

快速运行保留全部 36 个组合，但每个组合只做 **4 批、每批 512 条路径**；它用于检查流程，**不能复现完整研究的 Monte Carlo 结果表**。程序拒绝写入非空输出目录；再次运行时请选择新目录，以保留已有参照结果。[requirements-recorded.txt](requirements-recorded.txt) 记录原始环境；库版本和机器不同，可能影响随机数流、舍入及运行时间。

## 实验比较什么

三种状态方案都使用相同的梯形积分：

- **Exact**：通过非中心卡方分布模拟精确 CIR 网格转移。
- **Full truncation**：全截断 Euler；辅助状态允许为负，但漂移、扩散和折现使用其正部。
- **QE**：匹配条件均值和方差的二次—指数转移方案。

主实验使用三个合成参数情景、四种网格（每年 2、4、12、48 步）和五年期限。每个组合有 **80 个独立批次，每批 8,000 条路径**，合计 2,880 批、2,304 万条路径，种子基数为 `20261011`。年度网格下的价格抵消例子是单独的解析计算。

精确网格转移的仿射公式给出积分偏差参照。Monte Carlo 实验检查的是**价格及价格区间的覆盖率**；敏感度例子比较连续时间与离散积分的解析敏感度，未验证 Euler 或 QE 的 Monte Carlo Greek 估计器。

## 一组值得核查的结果

在合成的 `boundary` 情景和半年步长下：

| 方案 | 估计价格误差，单位本金 bp | 平均误差的 MCSE，bp | 连续时间价格覆盖次数 |
| --- | ---: | ---: | ---: |
| Exact | +1.397 | 1.114 | 76/80 |
| QE | +2.386 | 1.095 | 78/80 |
| Full truncation | −184.408 | 1.324 | 0/80 |

这一网格的解析积分偏差只有 **+1.219 bp**。全截断方案的单批区间平均半宽为 **23.546 bp**，仅报告采样误差条会遗漏主要的状态离散偏差。在固定的有偏网格上增加路径数不会消除该偏差。全部结果见 [summary.csv](summary.csv) 和[误差审计图](error_audit.png)。

覆盖率本身也有估计误差：若真实覆盖率为 95%，80 次重复的覆盖率 MCSE 约为 **2.44 个百分点**，不能用几个百分点的差异建立可靠排名。0/80 也不代表真实覆盖率严格为零；其 Wilson 95% 上限约为 4.58%。

## 单位与适用范围

- **价格误差 bp** 定义为单位本金下的 $10^4(\widehat P-P)$，**不是收益率基点**。
- **$G$** 是固定 $\kappa,\theta,\sigma$ 时对初始短率 $r_0$ 的敏感度，**不是整条收益率曲线平移的 DV01**。为保留参照数据，CSV 中的 `*_r0_dv01` 列名沿用原命名，存储的是 $10^{-4}G$，即初始短率上升一个基点时的一阶价格下降量。
- **平均误差的 MCSE** 描述独立批次间估计平均误差的不确定性，与单批价格置信区间不同。
- **覆盖率** 是“均值 ± 1.96 × 样本标准误”区间覆盖连续时间解析价格的比例，并非有限样本覆盖保证。

低 Feller 参数是边界压力情景，并未经过市场校准。方案按相同路径数比较，并非相同运行时间预算。参数估计不确定性、市场校准和模型设定错误不在本实验范围内。

## 文件导航

| 文件 | 内容 |
| --- | --- |
| [study.py](study.py) | 模拟、仿射参照、检查、结果表和绘图 |
| [test_validation.py](test_validation.py) | 自动化实现检查 |
| [report.md](report.md) / [report.zh-CN.pdf](report.zh-CN.pdf) | 推导、实验细节与局限 |
| [summary.csv](summary.csv) / [replicates.csv](replicates.csv) | 汇总结果／逐批次估计 |
| [analytic_grid.csv](analytic_grid.csv) | 精确网格定价、敏感度与积分偏差 |
| [cancellation_counterexamples.csv](cancellation_counterexamples.csv) | 价格与敏感度抵消反例 |
| [results.json](results.json) | 参数、种子、环境、检查与参照结果 |
| [requirements.txt](requirements.txt) / [requirements-recorded.txt](requirements-recorded.txt) | 安装依赖／原始环境版本记录 |

## 与已有工作的关系

离散积分递推是已有结果的透明复推，直接近邻是 **Peng 和 Kou（2008）第 2.5.2 节、Lemma 2**。离散化造成 Monte Carlo 区间覆盖失真，也已由 **Detemple、Garcia 和 Rindisbacher（2003）** 等研究。这个仓库将相关机制组织成可检查、可复现的定价与敏感度验证实验。

- Cox、Ingersoll 和 Ross（1985）：[A Theory of the Term Structure of Interest Rates](https://doi.org/10.2307/1911242)。
- Peng 和 Kou（2008）：[Connecting the Top-down to the Bottom-up: Pricing CDO under a Conditional Survival (CS) Model](https://english.phbs.pku.edu.cn/uploadfile/2019/0920/20190920023124792.pdf)。
- Detemple、Garcia 和 Rindisbacher（2003）：[Asymptotic Properties of Monte Carlo Estimators of Diffusion Processes](https://cirano.qc.ca/files/publications/2003s-11.pdf)。

数值方案的参考文献和后续研究方向见[技术报告](report.md)。

开发辅助：初稿及代码使用了 AI 辅助，仓库包含复现检查。
