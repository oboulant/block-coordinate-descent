# block-coordinate-descent
Implementation of Bleakley and Vert 2011

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

