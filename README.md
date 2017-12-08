# 1D-Lagrangian-hydrocode

This lagrangian code follows the adiabatic expansion of a hot spherical bubble in a uniform ambient
medium. The present setup is adjusted to spherical geometry.

This code is based on the Fortran version included in:
			Numerical methods in astrophysics: an Introduction - Peter Bodenheimer


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and usage purposes.

<!-- ### Prerequisites
 -->
<!-- What things you need to install the software and how to install them

```
Give examples
```-->


### Installing

First, clone the repository:

- Create a directory ```my_repo``` on your local machine 
- run: ```git clone https://github.com/eliotayache/1D-Lagrangian-hydrocode my-repo```

The main directory contains four sub-folders:

- C : Contains the lagrangian hydro source code in C
- C_Bondi : Same code with Bondi accretion feature
- Fortran : Contains the lagrangian hydro source code in Fortran
- Python : Contains plotting procedures

You can compile the code in the C and C_Bondi directories by typing:
```
cd C_Bondi 		# moving to C_bondi directory
gcc -o lh1 lh1.c
```

You now have an up and running version of the code!

## Using the code

Dynamical evolution is computed either in ```C``` or ```C_Bondi```. Plotting is done by using the
```plot.py``` python script located in the directory ```Python```.

### Compute dynamical evolution

Input values are set up in ```lh1.dat```. After updating these values, just run the code by typing:
```
./lh1
```

Output is saved in ```lh1.out``` and ```lh2.out```.

### Plotting results

Explain what these tests test and why

```
Give an example
```

<!-- ## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 
 -->

## Authors

* **Eliot Ayache** - *Initial work* - [github](https://github.com/eliotayache)

<!-- See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.
 -->

<!-- ## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
 -->
 
## Acknowledgments

This code is based on the Fortran version included in:
				Numerical methods in astrophysics: an Introduction - Peter Bodenheimer

