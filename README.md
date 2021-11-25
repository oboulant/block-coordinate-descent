# block-coordinate-descent

Implementation in C of the exact solution by Block Coordinate Descent for fast detection of multiple change-points.

It follows the following paper :

* J.-P. Vert and K. Bleakley, "Fast detection of multiple change-points shared by many signals using group LARS", In J. Lafferty, C. K. I. Williams, J. Shawe-Taylor, R.S. Zemel and A. Culotta (Eds), Advances in Neural Information Processing Systems 23 (NIPS), p.2343-2351, 2010. [[paper]](https://members.cbio.mines-paristech.fr/~jvert/svn/ngs/Lasso/article/groupLARS/nips2010/nips2010.pdf) [[supplementary informations]](https://members.cbio.mines-paristech.fr/~jvert/svn/ngs/Lasso/article/groupLARS/nips2010/supplementary.pdf) [[poster]](https://members.cbio.mines-paristech.fr/~jvert/publi/nips2010poster/poster.pdf)

The Matlab implementation by the authors can be found [here](https://members.cbio.mines-paristech.fr/~jvert/svn/GFLseg/html/).

Details about the pseudo code by the same authors can be found [here](https://hal.archives-ouvertes.fr/hal-00602121).

## Examples

### Clean and Compile

```
> make clean all
```

### Run

```
> ./main
```

### Change input data

There are two files in the folder `./data/` :

* `test1.txt`
* `test2.txt` 

If you want to change input data file, please 

* edit `./examples/main.c` with the proper file path,
* change `nb_lines` and `nb_cols` in accordance with you data.

## Tests

### Clean and Compile

```
> make clean test
```

### Run

```
> ./test
```

