# HurricaneSandy

**Microseismic Sources during Hurricane Sandy**

In this research, we studied the secondary microseismic signals excited by Hurricane Sandy recorded at USArray stations in the east of U.S. In one-hour time window, we found these signals were only correlated among stations aligned along close azimuths from the hurricane center. The travel-time differences were measured between the correlated stations, and these measurements were used to determine the source locations. We attributed these correlated seismic signals to two types of seismic sources, with one group of the seismic signals from the hurricane center and the other from coastal region. We further developed a hurricane seismic source model, to quantitatively describe the coupling among sea level pressure fluctuations, ocean waves and solid earth in the region of hurricane center. We also identified a strong seismic source near the coastal region in New England after Sandy’s dissipation, possibly related to subsequent storm surge in the area.

**Peer-reviewed Publication**
Chen, X., Tian, D., & Wen, L. (2015). Microseismic sources during Hurricane Sandy. Journal of Geophysical Research: Solid Earth, 120(9), 6386–6403. doi:10.1002/2015JB012282. 

**report of Science News Magazine**
https://www.sciencenews.org/article/hurricane%E2%80%99s-tiny-earthquakes-could-help-forecasters 

```src/cal_omp.c```: (OpenMP parallel computation) Source code of computing cross-correlation between two stations, and determining the travel-time differences of the seismic signals from seismic source to stations.

```src_locate/project.c```: Source code of locating source positions, using the results given by cal_omp

**项目描述**：

飓风是极具破坏力的自然现象之一，目前传统的飓风监测方法依赖于卫星、飞机与地面气象观测，而这些监测手段都具有很大的局限性。比如卫星不能很好地观测到飓风强度的急剧变化，且当飓风丧失一些热带气旋特征时，卫星观测手段就难以估计其强度。而气象监测部门的侦查飞机则无法对飓风内部的结构进行有效、持续和全面的观测，地面气象观测台站更是远离飓风中心，无法直接监测飓风中心的强度和结构。 

但是，飓风激发的海浪在海底会产生压强扰动，将能量传递给固体地球，产生微弱的地震波。通过提取地震信号，迅速定位震源位置并测量其强度，从而使有效、实时监测飓风成为可能。 

我们通过分析地震台站记录到的数据，研究建立了一个大气、海洋和固体地球相互耦合的物理模型，利用地震信号成功追踪2012年桑迪飓风中心位置，并监测飓风中心的海面气压扰动与海浪高度。研究还发现了传统卫星监测方法所无法观测到的、飓风消失以后的潜在灾害。研究显示美国新英格兰附近海域在桑迪飓风消失后仍然存在地震源，表明这些地区依然存在被海浪袭击的风险。 

这一研究成果表明，地震学可以成为一种有效、实时的飓风监测手段，为飓风内部动力学研究实时提供其中心区域气压扰动数据，同时监测飓风消失后依然存在的潜在灾害。该研究为现代地震学的应用提供了一个崭新的方向。

**项目职责**：

我是这项研究的主要负责人。我负责编写程序处理地震台站记录到的数据，分析并提出物理模型，撰写英文学术论文并发表。

**项目业绩**：

研究成果发表在英文SCI期刊： 
Chen, X., Tian, D., & Wen, L. (2015). Microseismic sources during Hurricane Sandy. Journal of Geophysical Research: Solid Earth, 120(9), 6386–6403. https://doi.org/10.1002/2015JB012282 

这项研究被国际Science News 专题新闻科普文章报道： 
https://www.sciencenews.org/article/hurricane%E2%80%99s-tiny-earthquakes-could-help-forecasters 

此外，国内各大媒体也对本研究成果进行了报道： 
http://seis.ustc.edu.cn/report/hurrican-sandy/ 
http://news.ustc.edu.cn/xwbl/201509/t20150925_227934.html 

项目的中文介绍页面： 
http://seis.ustc.edu.cn/research/hurricane-sandy/

**源码说明**

src/cal_omp.c:应用openmp并行计算不同台站对地震信号的互相关，并从中提取出地震信号到两个台站之间的走时差，并存储。
src_locate/project.c:使用cal_omp计算出的台站间地震波走时差，定位飓风激发的地震源位置
