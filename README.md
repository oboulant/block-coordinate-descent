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

## Run on your own data : Use the Python wrapper

If you want to use the present implementation on your own data, the easiest way is to use our Python wrapper. It can be found here : [https://github.com/oboulant/tamis](https://github.com/oboulant/tamis). 
It brings :

* The benefit of Python flexibility in terms of formating input data, 
* The benefit of C speed implementation. 

## Tests

We implemented unit tests using [Infer](https://fbinfer.com/docs/getting-started/). You can compile and run those tests using the following instructions. 

### Clean and Compile

```
> make clean test
```

### Run

```
> ./test
```

