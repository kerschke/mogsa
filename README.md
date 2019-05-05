# mogsa: Visualizing and Optimizing Multi-Objective Optimization Problems

## Description

Although many (real-world) optimization problems actually contain multiple objectives and
hence actually are multi-objective optimization problems, people often prefer to handle
them as single-objective problems (as those are easier to grasp).

In an attempt to make multi-objective problems more conceivable, *mogsa* provides the
means to visualize the landscapes of bi-objective problems -- more precisely: of its
gradients. As exemplarily shown by the heatmaps below, the visualizations reveal the
interactions between the local optima of the problem's different objectives in the
search and decision space.

![Examplary heatmap](https://raw.githubusercontent.com/kerschke/mogsa/master/images/heatmap_examples.png)

Based on the insights gained from these visualizations, we developed a multi-objective
gradient-based local search algorithm: the **m**ulti-**o**bjective **g**radient
**s**liding **a**lgorithm (MOGSA), which is also implemented within this package.


## Installation

Currently, *mogsa* is only available within this development version, however, we of
course also plan to submit it to CRAN in the near future.

In the mean time, feel free to use the development version of this package:

```r
install.packages("devtools")
devtools::install_github("kerschke/mogsa")
```


## Quickstart

To be continued ...


## Citation

If you like our proposed gradient-field heatmaps, cite our [EMO 2017 paper](http://link.springer.com/chapter/10.1007/978-3-319-54157-0_23) in your publications. Similarly, please refer to our [EMO 2019 paper](https://link.springer.com/chapter/10.1007/978-3-030-12598-1_11) when referencing MOGSA.

```
@inproceedings{KerschkeGrimme2017Expedition,
  author    = {Pascal Kerschke and Christian Grimme},
  title     = {{An Expedition to Multimodal Multi-Objective Optimization Landscapes}},
  booktitle = {{Proceedings of the 9$^{th}$ International Conference on Evolutionary Multi-Criterion Optimization (EMO)}},
  pages     = {329~--~343},
  series    = {{Lecture Notes in Computer Science (LNCS)}},
  volume    = {11411},
  editor    = {Heike Trautmann and Günter Rudolph and Kathrin Klamroth and Oliver Schütze and Margaret Wiecek and Yaochu Jin and Christian Grimme},
  year      = {2017},
  publisher = {Springer},
  address   = {M{\"u}nster, Germany},
  isbn      = {978-3-319-54157-0},
  doi       = {10.1007/978-3-319-54157-0_23},
  url       = {http://link.springer.com/chapter/10.1007/978-3-319-54157-0_23}
}

@inproceedings{GrimmeKerschkeTrautmann2019Multimodality,
  author    = {Christian Grimme and Pascal Kerschke and Heike Trautmann},
  title     = {{Multimodality in Multi-Objective Optimization --- More Boon than Bane?}},
  booktitle = {{Proceedings of the 10$^{th}$ International Conference on Evolutionary Multi-Criterion Optimization (EMO)}},
  pages     = {126~--~138},
  series    = {{Lecture Notes in Computer Science (LNCS)}},
  volume    = {11411},
  editor    = {Kalyanmoy Deb and Erik Goodman and Coello Carlos A. Coello and Kathrin Klamroth and Kaisa Miettinen and Sanaz Mostaghim and Patrick Reed},
  year      = {2019},
  publisher = {Springer},
  location  = {{East Lansing}, {MI}, {USA}},
  doi       = {10.1007/978-3-030-12598-1_11},
  url       = {https://link.springer.com/chapter/10.1007/978-3-030-12598-1_11}
}
```


## Contact

If you have any suggestions or ideas, or run into problems while running the code, please
use the [issue tracker](https://github.com/kerschke/mogsa/issues) or send me an e-mail <kerschke@uni-muenster.de>.


