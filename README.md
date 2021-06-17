# block-coordinate-descent
Implementation of Bleakley and Vert 2011

## Build

### Compile

```
> make all
```

### Clean workspace

```
> make clean
```

### Clean and compile

```
> make clean all
```

## Run

```
> ./main
```

## Change input data

There are two files in the folder `./data/` :

* `test1.txt`
* `test2.txt` 

If you want to change input data file, please 

* edit `./examples/main.c` with the proper file path,
* change `nb_lines` and `nb_cols` in accordance with you data.
