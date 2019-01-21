# 通过 Python API 使用 QIIME 2

## 什么是 QIIME2？

QIIME 2 是一款强大、可扩展和去中心化的微生物组分析包，强调数据分析透明。QIIME 2 可以使研究者从原始 DNA 序列开始分析，直接获取出版级的统计和图片结果。

### QIIME 2 原理图

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1fz8txisgcjj20kq0gltdj.jpg)

### QIIME 2 的各类界面

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1fz8txk4hvgj20u009t0yt.jpg)

## 准备工作

首先需要安装 QIIME2，直接用 conda 安装即可：

```
wget https://data.qiime2.org/distro/core/qiime2-2018.11-py35-linux-conda.yml
conda env create -n qiime2-2018.11 --file qiime2-2018.11-py35-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2018.11-py35-linux-conda.yml
```

激活环境：

```
source activate qiime2-2018.11
```

打开 Jupyter Notebook：

```
jupyter notebook
```

## 一些简单的示例

例如使用`rarefy`方法，输入文件为`q2-feature-table`插件产生的特征表，输出文件为样本深度一致的特征表。

```
>>> from qiime2.plugins import feature_table
>>> from qiime2 import Artifact
>>> unrarefied_table = Artifact.load('table.qza')
>>> rarefy_result = feature_table.methods.rarefy(table=unrarefied_table, sampling_depth=100)
>>> rarefied_table = rarefy_result.rarefied_table
```

比起`cli`的方便之处就在于，这样就可以直接作为 `Python` 对象来访问基础数据。例如，创建为`biom.Table`对象的`rarefied_table`：

```
>>> import biom
>>> biom_table = rarefied_table.view(biom.Table)
>>> print(biom_table.head())
# Constructed from biom file
#OTU ID      L1S105  L1S140  L1S208  L1S257  L1S281
b32621bcd86cb99e846d8f6fee7c9ab8     25.0    31.0    27.0    29.0    23.0
99647b51f775c8ddde8ed36a7d60dbcd     0.0     0.0     0.0     0.0     0.0
d599ebe277afb0dfd4ad3c2176afc50e     0.0     0.0     0.0     0.0     0.0
51121722488d0c3da1388d1b117cd239     0.0     0.0     0.0     0.0     0.0
1016319c25196d73bdb3096d86a9df2f     11.0    17.0    12.0    4.0     2.0
```

除此之外，也可以直接将他们看做`Pandas`对象：

```
>>> import pandas as pd
>>> df = rarefied_table.view(pd.DataFrame)
>>> df.head()
        b32621bcd86cb99e846d8f6fee7c9ab8  99647b51f775c8ddde8ed36a7d60dbcd  \
L1S105                              25.0                               0.0
L1S140                              31.0                               0.0
L1S208                              27.0                               0.0
L1S257                              29.0                               0.0
L1S281                              23.0                               0.0
...
```

QIIME 的一个强大功能就是，你既可以从 QIIME2 中导出不同的类型来进行处理，也可以将这些数据处理后导入QIIME（但这也有个缺点，QIIME会无法跟踪你的数据来源）。比如，你可以直接导入为`Pandas`对象来处理。

```
imported_artifact = Artifact.import_data("FeatureTable[Frequency]", df)
```

`rarefied_table` 可以直接传递给其他 QIIME 2 插件。使用 `q2-diversity` 插件计算 alpha 多样性`observed_otus`（OTU数量）指数。生成的对象将是`SampleData [AlphaDiversity]`类型，同样可以用`pd.Series`来进行处理。

```
>>> from qiime2.plugins import diversity
>>> alpha_result = diversity.methods.alpha(table=rarefied_table, metric='observed_otus')
>>> alpha_diversity = alpha_result.alpha_diversity
>>> alpha_diversity.view(pd.Series)
L1S105    24
L1S140    19
L1S208    25
L1S257    30
L1S281    29
L1S57     23
L1S76     20
L1S8      17
...
Name: observed_otus, dtype: int64
```

最后，我们可以将这些对象保存为 QIIME 可以识别的`.qza`文件，并退出：

```
>>> rarefied_table.save('rare.qza')
'rare.qza'
>>> alpha_diversity.save('oo.qza')
'oo.qza'
>>> exit
```

## Moving Pictures tutorial

下面，就开始用 Python 完整地走一遍 QIIME 2 的入门教程吧~

在本教程中，你将使用 QIIME 2 在五个时间点对来自四个身体部位的两个人的微生物组样本进行分析，第一个时间点紧接着是抗生素的使用。基于这些样本的研究文章《Moving pictures of the human microbiome》在2011年发表于Genome Biology。本教程中使用的数据基于Illumina HiSeq产出，使用地球微生物组项目扩增16S rRNA基因高变区4（V4）测序的方法。

### 准备工作

```
# 创建本节学习目录
mkdir qiime2-moving-pictures-tutorial
cd qiime2-moving-pictures-tutorial
# 下载 metadata ，需要翻墙
wget \
  -O "sample-metadata.tsv" \
  "https://data.qiime2.org/2018.11/tutorials/moving-pictures/sample_metadata.tsv"
# 下载测序数据
mkdir -p emp-single-end-sequences
wget -O "emp-single-end-sequences/barcodes.fastq.gz" "https://data.qiime2.org/2018.11/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz"
wget -O "emp-single-end-sequences/sequences.fastq.gz" "https://data.qiime2.org/2018.11/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz"
# 打开 jupyter notebook
jupyter notebook
```

- 导入模块

```python
%matplotlib inline
import qiime2
from tempfile import mkdtemp
from qiime2.plugins import demux, deblur, quality_filter, \
                           metadata, feature_table, alignment, \
                           phylogeny, diversity, emperor, feature_classifier, \
                           taxa, composition
```

- 导入测序数据

首先我们需要导入数据，这个对象的语义类型是`EMPSingleEndSequences`。`emp-single-end-sequences`文件夹中包括`barcodes.fastq.gz`和 `sequences.fastq.gz`，`sequences.fastq.gz`是包含多样本混合的序列文件，所以需要`barcode`来拆分样本。

```python
single_end_sequences = qiime2.Artifact.import_data('EMPSingleEndSequences', './emp-single-end-sequences')
```

- 导入 metadata

```python
sample_metadata = qiime2.Metadata.load('./sample-metadata.tsv')
```

### 拆分样品

在`sample-metadata.tsv`中包含了每个样本对应的`barcode`，你可以运行以下命令来对序列进行样本拆分。（`demux emp-single`命令指的是这些序列是根据地球微生物学方法添加的条形码，并且是单端序列）。`demux_sequences`对象包含样本拆分后的序列。


```python
demux_sequences = demux.methods.emp_single(single_end_sequences,
                                           sample_metadata.get_column('BarcodeSequence'))
```

在样本拆分之后，你可以生成拆分结果的摘要，可以看到每个样本有多少序列，还可以获得序列数据中每个位置处序列质量分布的摘要。

```python
demux_summary = demux.visualizers.summarize(demux_sequences.per_sample_sequences)
demux_summary.visualization.export_data('./01-demux_summary')
```

这里我用`demux_summary.visualization.export_data`来保存网页（QIIME2的结果是以网页来展示的）到本地，打开文件夹中的`index.html`即可查看结果。其实有好几种方式可以显示结果：

1. 导出`.qzv`用`cli`来看结果：

```
demux_summary.visualization.save('demux_summary.qzv')
!qiime tools view demux_summary.qzv
```

2. 在`jupyter notebook`中简单粗暴地插入网页：

```
from IPython.core.display import  HTML
HTML('<iframe height = 1000, width = 1000, src = "{}/index.html"> </iframe>'\
    .format('./01-demux_summary'))
```

### 序列质控和生成特征表

QIIME 2 插件多种质量控制方法可选，包括 DADA2 、Deblur 和基于基本质量分数的过滤。在本教程中，我们使用 DADA2 和 Deblur 介绍这个步骤。这些步骤是可互相替换的，因此你可以使用自己喜欢的方法。这两种方法的结果会产生两个对象 `FeatureTable[Frequency]`和`FeatureData[Sequence]`，`Frequency`对象包含数据集中每个样本中每个唯一序列的计数（频率），`Sequence`对象将`FeatureTable`中的特征ID与序列对应。

#### Option 1: DADA2

##### 质控

```python
#skipped
from qiime2.plugins import dada2
qual_control = dada2.actions.denoise_single(demultiplexed_seqs=demux_sequences.per_sample_sequences, trunc_len=120, trim_left=0)
```

`dada2 denoise-single`方法需要两个用于质量过滤的参数：`--p-trim-left m`，它去除每个序列的前m个碱基（如引物、barcode）；`--p-trunc-len n`，它在位置n截断每个序列。这允许用户去除序列的低质量区域、引物等。需要查看上面`demux_summary`的结果图以确定这两个值。

##### 对特征表统计进行进行可视化

```python
stats_meta = metadata.actions.tabulate(input=qual_control.denoising_stats.view(qiime2.Metadata))
stats_meta.visualization.export_data('./02-stats_meta')
ft_sum = feature_table.actions.summarize(table = qual_control.table)
ft_sum.visualization.export_data('./03-ft_sum')
tab_seqs = feature_table.actions.tabulate_seqs(data=qual_control.representative_sequences)
tab_seqs.visualization.export_data('./04-tab_seqs')
```

#### Option 2: Deblur

Deblur使用序列错误配置文件将错误的序列与从其来源的真实生物序列相关联，从而得到高质量的序列数据，主要为两个步骤。首先，使用基于质量分数的初始质量过滤：

```python
demux_filter_stats = quality_filter.methods.q_score(demux_sequences.per_sample_sequences)
```

接下来，使用`deblur denoise-16S`方法。此方法需要一个用于质量过滤的参数，即截断位置n长度的序列的`trim_length`。通常，Deblur开发人员建议将该值设置为质量分数中位数开始下降至低质量区时的长度。在本次数据上，质量图表明合理的选择应该是在 115 至 130 序列位置范围内。这里我们也使用`trim_length=120`参数。

```python
deblur_sequences = deblur.methods.denoise_16S(demux_sequences.per_sample_sequences,
                                              trim_length=120,
                                              sample_stats=True)
filter_stats = metadata.visualizers.tabulate(demux_filter_stats.filter_stats.view(qiime2.Metadata))
filter_stats.visualization.export_data('./02-filter_stats')
deblur_stats = deblur.visualizers.visualize_stats(deblur_sequences.stats)
deblur_stats.visualization.export_data('./03-deblur_stats')
```

### 特征表和特征序列汇总

在质控完成之后，可以使用以下两个命令来探索数据结果，这两个命令将创建数据的可视化摘要。特性表汇总命令（`feature_table.visualizers.summarize`）将向你提供关于与每个样品和每个特性相关联的序列数量、这些分布的直方图以及一些相关的汇总统计数据的信息。特征表序列表格`feature_table.visualizers.tabulate_seqs`命令将提供特征 ID 到序列的映射，并提供链接以针对 NCBI nt 数据库 BLAST 每个序列。（下面只用 Deblur 方法的结果来示例）

```python
output_viz = feature_table.visualizers.summarize(deblur_sequences.table)
output_viz.visualization.export_data('./04-output_viz')
deblur_feature_table_summary = feature_table.visualizers.tabulate_seqs(deblur_sequences.representative_sequences)
deblur_feature_table_summary.visualization.export_data('./05-deblur_feature_table_summary')
```

### 构建进化树用于多样性分析

QIIME 2支持几种系统发育多样性度量方法，包括`Faith’s Phylogenetic Diversity`、`weighted`和`unweighted UniFrac`。除了每个样本的特征计数（即 QIIME2 对象`FeatureTable[Frequency]`）之外，这些度量还需要将特征彼此关联结合有根进化树。此信息将存储在一个 QIIME 2 对象的有根系统发育对象`Phylogeny[Rooted]`中。为了生成系统发育树，我们将使用`q2-phylogeny`插件中的`align-to-tree-mafft-fasttree`工作流程。

首先，工作流程使用 mafft 程序执行对`FeatureData[Sequence]`中的序列进行多序列比对，以创建 QIIME 2对象`FeatureData[AlignedSequence]`。

```python
mafft_alignment = alignment.methods.mafft(deblur_sequences.representative_sequences)
```

接下来，流程屏蔽（mask 或过滤）对齐的的高度可变区(高变区)。这些位置通常被认为会增加系统发育树的噪声。

```python
masked_mafft_alignment = alignment.methods.mask(mafft_alignment.alignment)
```

随后，流程应用 FastTree 基于过滤后的比对结果生成系统发育树。

```python
unrooted_tree = phylogeny.methods.fasttree(masked_mafft_alignment.masked_alignment)
```

FastTree 程序创建的是一个无根树，因此在本节的最后一步中，应用根中点法将树的根放置在无根树中最长端到端距离的中点，从而形成有根树。

```python
rooted_tree = phylogeny.methods.midpoint_root(unrooted_tree.tree)
```

### Alpha和beta多样性分析

**Alpha and beta diversity analysis**

QIIME 2 的多样性分析可以通过`diversity`插件实现，该插件支持计算 α 和 β 多样性指数、并应用相关的统计检验以及生成交互式可视化图表。

首先使用`core_metrics_phylogenetic`方法，该方法将`FeatureTable[Frequency]`抽平到用户指定的测序深度，计算常用几种 α 和 β 多样性指数，并使用`Emperor`为每个 β 多样性指数生成主坐标分析（PCoA）图。默认情况下计算的方法有：

- α多样性
    - 香农（Shannon's）多样性指数（群落丰富度的定量度量，即包括丰富度`richness`和均匀度`evenness`两个层面）
    - Observed OTUs（群落丰富度的定性度量，只包括丰富度）
    - Faith’s系统发育多样性（包含特征之间的系统发育关系的群落丰富度的定性度量）
    - 均匀度（或 Pielou’s均匀度；群落均匀度的度量）
- β多样性
    - Jaccard距离（群落差异的定性度量，即只考虑种类，不考虑丰度）
    - Bray-Curtis距离（群落差异的定量度量）
    - 非加权UniFrac距离（包含特征之间的系统发育关系的群落差异定性度量）
    - 加权UniFrac距离（包含特征之间的系统发育关系的群落差异定量度量）

需要提供的一个重要参数是`sampling_depth`，它是指定重采样（即稀疏/稀释）深度。因为大多数多样指数对不同样本的不同测序深度敏感，所以这个脚本将随机地将每个样本的测序量重新采样至该参数提供的值。例如，如果提供`sampling_depth 500`，则此步骤将对每个样本中的计数进行无放回抽样，从而使得结果表中的每个样本的总计数为 500 。如果任何样本的总计数小于该值，那么这些样本将从多样性分析中删除。建议通过查看上面的`output_viz`对象中所呈现的信息并**选择一个尽可能高的值（因此每个样本保留更多的序列）同时尽可能少地排除样本来进行选择**。

```python
core_metrics = diversity.pipelines.core_metrics_phylogenetic(table = deblur_sequences.table,
                                                             phylogeny = rooted_tree.rooted_tree,
                                                             sampling_depth = 1109,
                                                             metadata = sample_metadata)
```

这里，我们将`sampling_depth`参数设置为 1109 。这个值是根据`L3S341`样本中的序列数量来选择的，因为它与接下来几个序列计数较高的样本中的序列数量接近，并且因为它仅比序列较少的一个样本中的序列数量高。这将允许我们保留大部分样品。具有较少序列的一个样本将从`core-metrics-phylogenetic`分析和任何使用这些结果的下游分析中删除。

### Alpha多样性组间显著性分析和可视化

在计算多样性之后，我们可以开始探索样本的微生物组成。首先将`metadata`和 alpha 多样性数据相关联。下面以 faith-pd 为例探索不同元数据条件下的组间差异：

```python
faith_pd_group_significance = diversity.actions.alpha_group_significance(core_metrics.faith_pd_vector,
                                                                         sample_metadata)
faith_pd_group_significance.visualization.export_data('./06-faith_pd_group_significance')
```

```python
evenness_group_significance = diversity.actions.alpha_group_significance(core_metrics.evenness_vector,
                                                                         sample_metadata)
evenness_group_significance.visualization.export_data('./07-evenness_group_significance')
```

### beta group significance

接下来，我们将使用 PERMANOVA 方法 beta-group-significance 分析分类元数据背景下的样本组合。以下命令将测试一组样本之间的距离，是否比来自其他组（例如，舌头、左手掌和右手掌）的样本彼此更相似，例如来自同一身体部位（例如肠）的样本。这里，使用两个示例元数据列将此应用到未加权的 UniFrac 距离，如下所示：

```python
uUniFrac_BodySite_significance = diversity.actions.beta_group_significance(core_metrics.unweighted_unifrac_distance_matrix,
                                                                           sample_metadata.get_column('BodySite'))
uUniFrac_BodySite_significance.visualization.export_data('./08-uUniFrac_BodySite_significance')
```

```python
uUniFrac_Subject_significance = diversity.actions.beta_group_significance(core_metrics.unweighted_unifrac_distance_matrix,
                                                                          sample_metadata.get_column('Subject'))
uUniFrac_Subject_significance.visualization.export_data('./09-uUniFrac_Subject_significance')
```

最后，分类是在样本元数据分组间探索微生物群落组成差异的流行方法。可以使用 Emperor 工具在示例元数据背景下探索主坐标（PCoA）绘图。虽然`core_metrics_phylogenetic`命令已经生成了一些 Emperor 图，但还可以传递一个可选的参数`custom_axes`来探索时间序列数据。我们将为未加权的 UniFrac 的PCoA结果生成 Emperor 图，所得到的图将包含主坐标1、主坐标2和实验开始以来的天数（days since the experiment start）的轴。使用最后一个轴来探索这些样本是如何随时间变化的。


```python
emperor_plot = emperor.visualizers.plot(core_metrics.unweighted_unifrac_pcoa_results,
                                        sample_metadata,
                                        custom_axes=['DaysSinceExperimentStart'])
emperor_plot.visualization.export_data('./10-emperor_plot')
```

### Alpha 稀释曲线

在本节中，使用`diversity.actions.alpha_rarefaction`来探索 α 多样性与采样深度的关系。它在多个采样深度处计算一个或多个α多样性指数，范围介于1（可选`min_depth`控制）和作为`max_depth`提供值之间。在每个采样深度，将生成 10 个抽样表，并对表中的所有样本计算 alpha 多样性指数计算。

```python
rarefaction = diversity.actions.alpha_rarefaction(table = deblur_sequences.table,
                                                  max_depth = 4000,
                                                  phylogeny = rooted_tree.rooted_tree,
                                                  metadata = sample_metadata)
rarefaction.visualization.export_data('./11-rarefaction')
```

### 物种组成分析

在这一节中，我们将探索样本的物种组成，并将其与样本元数据再次组合。这个过程的第一步是为`FeatureData[Sequence]`的序列进行物种注释。我们将使用经过 Naive Bayes 分类器预训练的，并由`q2-feature-classifier`插件来完成这项工作。这个分类器是在`Greengenes 13_8 99% OTU`上训练的，其中序列被修剪到仅包括来自 16S 区域的 250 个碱基，该 16S 区域在该分析中采用 V4 区域的 515F/806R 引物扩增并测序。我们将把这个分类器应用到序列中，并且可以生成从序列到物种注释结果关联的可视化。

- 下载物种注释数据库制作的分类器：

```python
!mkdir ./classifier
!wget -O "./classifier/gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2018.11/common/gg-13-8-99-515-806-nb-classifier.qza"
```

- 导入分类器：

```python
gg_classifier = qiime2.Artifact.import_data('TaxonomicClassifier', './classifier/')
```

- 分析物种组成：

```python
taxonomy = feature_classifier.methods.classify_sklearn(reads = deblur_sequences.representative_sequences,
                                                       classifier = gg_classifier)
taxonomy_classification = metadata.visualizers.tabulate(taxonomy.classification.view(qiime2.Metadata))
taxonomy_classification.visualization.export_data('./12-taxonomy_classification')
```

接下来，我们可以用交互式条形图查看样本的分类组成。使用以下命令绘图成堆叠柱状图，然后打开可视化。

```python
taxa_bar_plot = taxa.visualizers.barplot(deblur_sequences.table, taxonomy.classification, sample_metadata)
taxa_bar_plot.visualization.export_data('./13-taxa_bar_plot')
```

### 使用 ANCOM 进行差异丰度分析

ANCOM 可用于识别不同样本组中丰度差异的特征。ANCOM 是在`composition`插件中实现的。ANCOM 假设很少（小于约25%）的特征在组之间改变。**如果你期望在组之间有更多的特性正在改变，那么就不应该使用ANCOM**，因为它更容易出错（I类和II类错误都有可能增加）。因为我们预期身体部位的许多特征都会发生变化，所以在本教程中，我们将过滤完整的特征表后只包含肠道样本。然后，我们将应用 ANCOM 来确定哪种（如果有的话）序列变体在我们两个受试者的肠道样本中丰度存在差异。

- 首先创建一个只包含肠道样本的特征表：

```python
gut_deblur = feature_table.methods.filter_samples(deblur_sequences.table,
                                                  metadata = sample_metadata,
                                                  where = "BodySite='gut'")
```

ANCOM 基于每个样本的特征频率对`FeatureTable[Composition]`进行操作，但是不能容忍零。为了构建组成对象，必须提供一个添加伪计数（一种遗失值插补方法）的`FeatureTable[Frequency]`对象，这将产生`FeatureTable[Composition]`对象。

```python
gut_deblur_composition = composition.actions.add_pseudocount(gut_deblur.filtered_table)

ancom_gut_deblur = composition.actions.ancom(table = gut_deblur_composition.composition_table,
                                             metadata = sample_metadata.get_column('Subject'))
ancom_gut_deblur.visualization.export_data('./14-ancom_gut_deblur')
```

我们也经常对在特定的分类学层次上执行差异丰度检验。为此，我们可以在感兴趣的分类级别上折叠 FeatureTable[Frequency] 中的特性，然后重新运行上述步骤。下面，我们将特征表折叠到属级别（即 Greengenes 分类法的第6级）：


```python
gut_table_l6 = taxa.methods.collapse(table = gut_deblur.filtered_table,
                                     taxonomy = taxonomy.classification,
                                     level = 6)
gut_table_l6_composition = composition.actions.add_pseudocount(gut_table_l6.collapsed_table)

ancom_gut_table_l6 = composition.actions.ancom(table = gut_table_l6_composition.composition_table,
                                               metadata = sample_metadata.get_column('Subject'))
ancom_gut_table_l6.visualization.export_data('./15-ancom_gut_table_l6')
```

## Reference

- https://docs.qiime2.org/2018.11/tutorials/moving-pictures/
- https://docs.qiime2.org/2018.11/interfaces/artifact-api/
- https://gist.github.com/tkosciol/29de5198a4be81559a075756c2490fde
- https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition/tree/master/Chapter10
- https://forum.qiime2.org/t/python-api-visualizations/3097
- 宏基因组 - QIIME 2用户文档. 4人体各部位微生物组分析实战Moving Pictures(2018.11)


>PS：七牛云竟然挂了（是不是村网通了QAQ），找了个微博图床还蛮好用的。Chrome 商店搜索“新浪微博图床”即可~
