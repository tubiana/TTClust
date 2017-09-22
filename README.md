# TTClust : A molecular simulation trajectory clustering programm
---
## DESCRIPTION
TTclust is a python program used to clusterized molecular dynamics simulation trajectories. It requires just a trajectory and a topology file (compatible with most molecular dynamic packages such as Amber, Gromacs, Chramm, Namd or trajectory in PDB format thanks to the MDtraj package).
Easy to use, the program produces a visual feedback of the clustering through a dendrogram graph. Other graphic representations are made to describe all clusters (see OUTPUT EXAMPLE part for more details).

### Python Compatibility 
This program is compatible with python 2.7.x and 3.x.
You can have a warning for matplotlib with python 3.x but the program still works

### Dependancies and installation: 
Following packages are needed: 
  - argparse
  - argcomplete (for autocompletion, optional)
  - cython (for mdtraj)
  - mdtraj (version >= 0.17)
  - progressbar
  - datetime *(present in default python library)*
  - glob *(present in default python library)*
  - matplotlib
  - scipy (version >= 0.18)
  - prettytable
  - sklearn (version >= 0.18)
  
You will find a file **requirements.txt**. You can install all requiered 
package with this PIP command:  `sudo pip install -r requirements.txt`
Note : sometimes mdtraj fails to install. Please install mannualy cython before in this case `sudo pip install cython` then `sudo pip install -r requirements.txt`

#### For Mac user
If you have issues with pip, try first to add to pip the `--ignore-installed` argument : `sudo pip install --ignore-installed -r requirements.txt`
If it still doesn't work, it's maybe because of the System Integrity Protection (SIP).
I suggest you in this case to install ANACONDA or MINICONDA and restart your terminal afterwards. 
Normally, the pip command should work because your default python will be the anaconda (or miniconda) python

To activate autocompletion for the argpase module you have to use this command 
(only once) `sudo activate-global-python-argcomplete`

#### Atoms selection
For Selection syntax, use the one from MDTraj (http://mdtraj.org/latest/atom_selection.html).
You can specify different selections for the calculation: 
 - **st** is used to extract a part of the trajectory (if this one is too big).
 ..* leave this blank if you want to keep your trajectory intact.
 - **sa** used to align your trajectory on the same reference (eg: a chain or 
 ..* the backbone) 
 - **sr** is used for the clustering (rmsd calculated on the atom selected with 
 ..* this string)
 
#### NOTE on Nucleic Acids
MDTRAJ doesn't have nucleic acid keywords yet. We've implemented some keywords that will be altered to match DNA/RNA....
Keywords added : 
 - **dna** : selection based on the residue name (DA/DT/DC/DG)
 - **rna** : selection based on the residue name (A/T/G/C or RA/RT/RG/RC)
 - **backbone_na** : backbone of nucleic acid. Selection based on the residue name and atom name (P, O3', O5', C3', C4', C5')
 - **base** : selection base on the residue name and atom name. select RNA or DNA and exclude backbone_na, sugar atoms and hydrogen
 - **base_rna** : same as *base* but for RNA
 - **base_dna** : same as *base* but for DNA

 Theses selection keywords can be used with other MDTRAJ selection keywords, e.g.:
 - "protein and not dna"
 - "rna and not type H"

 
#### Clustering Methods
With the scipy module, several methods for clustering are available. Change the 
method used with the *-m* argument. Methods available are: 
 - single
 - complete
 - average
 - weighted
 - centroid
 - median
 - **ward** (DEFAULT)

3 possibilities  are available for the calculation: 

1. give the number of clusters you want. Eg: if you want 3 clusters, use the argument
..* **-ng 3**
2. give a cutoff for the clustering. The final clustering are made from a
..* dendrogram and this cutoff is used for the distance cutoff. If you want to
..* set this cutoff by hand, use the argument **-cc X** (X is the cutoff)
3. Choose your cutoff by clicking on the matplotlib windows (on the dendrogram)
..* in this case don't use the other arguments. **recommended for the first 
 clustering**

#### Distance Matrix
The distance matrix can be long to calculate depending on your trajectory size.
That's why this matrix is saved on the ".npy" format, in order to be used later.
The name of the matrix will be the name of your selection string for clustering (*sr*)
If you use the same selection string for clustering (*sr*) the matrix will be detected
and the programme will ask you if you want to use it again (Y), to recalculate this
matrix (N) or choose another matrix (O). If you want to use the saved matrix without
interactive this interactive question) add in argument **-i n** which will deactivate
the interactive prompt.



## ARGUMENTS 
```text
  -h, --help            show this help message and exit
  -f TRAJ, --traj TRAJ  trajectory file
  -t TOP, --top TOP     topfile
  -o OUTPUT, --output OUTPUT (default: clustering.log)
                        logfile 
  -st SELECT_TRAJ, --select_traj SELECT_TRAJ (default: all)
                        selection syntax for trajectory extraction, with QUOTE 
  -sa SELECT_ALIGNEMENT, --select_alignement SELECT_ALIGNEMENT (default: backbone)
                        selection syntax for alignement with QUOTE
						If you don't want alignement : use "none"
  -sr SELECT_RMSD, --select_rmsd SELECT_RMSD (default: backbone)
                        selection syntax for RMSD with QUOTE 
  -m METHOD, --method METHOD (default: ward)
                        method for clustering: single; complete; average;
                        weighted; centroid; median and ward
  -cc CUTOFF, --cutoff CUTOFF
                        cutoff for clusterization from hierarchical clusturing
                        with Scipy. If you choose to clic on the graph, cutoff
                        will be the clicked value 
  -ng NGROUP, --ngroup NGROUP
                        number of group wanted. Use the maxclust method to
                        clusterize in this case
  -i INTERACTIVE, --interactive INTERACTIVE
                        Interactive mode for distance matrix (Y/n)
```
    
## USAGE : 
There is some example usage with the examples files givent on the "example" folder. 
Please note that the trajectory is reduce to the backbone to reduce size of the git archive.
Caution : You have to put quote beside your selection string (for *sr*, *st*, and *sa* arguments) 
 - Simple usage (clustering on backbone, logfile is called clustering.log, output folder is "clustering")
 ```python ttclust.py -f examples/example.xtc -t examples/example.pdb```
 - Clustering on residues 30 to 200 and backbone
 ```python ttclust.py -f examples/example.xtc -t examples/example.pdb -sr "residue 30 to 200 and backbone" -l res30-200.log```
 - Clustering on CA atom and save this part of the trajectory with a cutoff of 2.75
 ```python ttclust.py -f examples/example.xtc -t examples/example.pdb -sr "name CA" -st "name CA" -cc 2.75 -l CA-c2.75.log```
 - Clustering on backbone of the protein and chain A (note that with mdtraj there is no chaine name, but chaine ID starting from 0) with 10 clusters only
 ```python ttclust.py -f examples/example.xtc -t examples/example.pdb -sr "protein and backbone and chainid 0" -l backbone-chainA.log -ng 10 ```
- Note For PDB trajectory, don't use the **-t** argument
```python TrajectoryClustering.py -f traj.pdb -st "protein" -sr "backbone"```


## OUTPUT EXAMPLE
#### Cluster result (structure)
PDBs are saved in a new folder Cluster_PDB.  
You can find in the PDB name the *cluster number*, *size* and the *representative frame*
of the cluster (*ie* the frame number of the saved structure)  
Example: C1-f11-s44.pdb corresponds to the cluter 1 made of 44 structures and the
saved frame (representative) is the frame 11.

#### Logfile
In the log file you will find all arguments given to the program with
cluster information:
 - **size**: number structure in the cluster
 - **representative frame**: frame with lowest RMSD between all other frame of the cluster
 - **Members**: all frame belonging to the cluster
 - **spread**: mean RMSD between all frame in the cluster
 - **RMSD between clusters**: A tab with the RMSD between clusters
 - **Average RMSD between clusters**: the average RMSD between clusters.

#### Dendrogram
A dendrogram is generated at the end of the clustering with the corresponding cluster colors.
The name of this file will be the same as the logfile with a ".png" extension 
example: example.log --> example.png
![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/backbone/backbone-den.png "Dendrogram example")  
The grey horizontal line is the cutoff value used.


#### LinearProjection representation
A linear projection of cluster is made for the trajectory.
![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/backbone/backbone-linear.png "linear-proj example")
Every barline represents a frame and the color a cluster number.
Note that : 
 - If less or equal than 12 clusters : a defined color map was made in this order :
   red, blue, lime, gold, darkorchid, orange, deepskyblue, brown, gray, black, darkgreen, navy
 - Else, the matplotlib "hsv" color map is used but the color change according to
   the number of clusters.

#### Barplot representation
A vertical barplot is generated to have a overview of the cluster size. the barcolor corresponds to the clusters color in the LinearProjection representation and dendrogram cluster's color.
![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/backbone/backbone-hist.png "histogram example")

#### 2D distance projection
A 2D projection of the distance(RMSD) between the representative frame of each cluster is made. The method used is the multimentional scaling method from the skilearn python module.
![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/backbone/backbone-dist.png "2D Distance example")
We can follow the evolution of each cluster thanks to the relative distance between them. The color of points is the same as for other graphs (ie. cluster colors) and the radius of each point depend on the cluster spread.

#### Distance matrix plot
A plot of the distance matrix is also made and allow to visualize the distance between two frames easily. 

![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/backbone/backbone-distmat.png "Distance Matrix plot example")


## Licence
This program is under the GNU GPLv3 licence, which means that anyone who 
distributes your code or a derivative work to make the source available under 
the same terms, and also provides an express grant of patent rights from 
contributors to users.
