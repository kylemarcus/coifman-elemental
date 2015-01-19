# Coifman Elemental

###Step 1:

[Install libElemental]

**Note:** Before install, load the following modules
```sh
$ moudle load cmake gcc
```

**Note:** Use the following cmake command (change install dir)
```sh
$ cmake -D CMAKE_INSTALL_PREFIX=/user/kmarcus2/Elemental/install -D CMAKE_CXX_COMPILER=/util/academic/gcc/gcc-4.8.2/bin/g++ -D CMAKE_C_COMPILER=/util/academic/gcc/gcc-4.8.2/bin/gcc ..
```

###Step 2:

You might have to change the some paths in the make file, then just run make
```sh
$ make
```

###Step 3:

Run the program (usually on a debug node)

```sh
$ mpirun -np 1 ./kyle-coifman-elemental Topsar360m.txt
```

###Config:

The configuraion file is **coifman.h**, change this and re-compile

[Install libElemental]:http://libelemental.org/documentation/dev/build.html#building-elemental
