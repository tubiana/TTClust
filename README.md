# TrajectorClustering
---
## DESCRIPTION
It uses the MDTraj package so you can use a lot of different type of files (
This program was developped in order to clusterize molecular dynamic trajectories.
Amber, gromacs, chamm, namd, PDB...)

### Python Compatibility 
This program is compatible with python 2.7.x and 3.x.
You can have a warning for matplotlib with python 3.x but the program still works

### Dependancies: 
Following packages are needed: 
  - argparse
  - argcomplete (for autocompletion, optional)
  - cython (for mdtraj)
  - mdtraj
  - progressbar
  - datetime *(present in default python library)*
  - glob *(present in default python library)*
  - matplotlib
  - scipy
  
You will find a file **requirements.txt**. You can install all requiered 
package wich this PIP:  `sudo pip install -r requirements.txt`

To activate autocompletion for the argpase module you have to use this command 
(only once) `sudo activate-global-python-argcomplete`

#### Atoms selection
For Selection syntax, use the one from MDTraj (http://mdtraj.org/latest/atom_selection.html).
You can specify different selections for the calculation: 
 - **st** is used to extract a part of the trajectory (if this one is too big).
 ..* leave this blank if you want to keep your trajectory intact.
 - **sa** used to align your trajectoriy on the same reference (eg: a chain or 
 ..* the backbone) 
 - **sr** is used for the clustering (rmsd calculated on the atom selected with 
 ..* this string)

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

3 possibilites are available for the calculation: 

1. give the number of clusters you want. Eg: if you want 3 clusters, use the argument
..* **-ng 3**
2. give a cutoff for the clustering. The final clustering are made from a
..* dendrogram and this cutoff is used for the distance cutoff. If you want to
..* set this cutoff by hand, use the argument **-cc X** (X is the cutoff)
3. Choose your cutoff by clicking on the matplotlib windows (on the dendrogram)
..* in this case don't use the other arguments. **recommended for the first 
 clustering**

#### Distance Matrix
The distance matrix can be long to calculated depending on your trajectory size.
That's why this matrix is saved on the ".npy" format, in order to be used later.
The name of the matrix will be the name of your selection string for clustering (*sr*)
If you use the same selection string for clustering (*sr*) the matrix will be detected
and the programe will ask you if you want to use it again (Y), to recalculate this
matrix (N) or choose an other matrix (O). If you want to use the saved matrix without
interactive this interactive question) add in argument **-i n** which will desactivate
the interactive prompt.



## ARGUMENTS 
```text
  -h, --help            show this help message and exit
  -f TRAJ, --traj TRAJ  trajectory file
  -t TOP, --top TOP     topfile
  -o OUTPUT, --output OUTPUT
                        logfile (logfile.txt)
  -st SELECT_TRAJ, --select_traj SELECT_TRAJ
                        selection syntax for trajectory extraction (all)
  -sa SELECT_ALIGNEMENT, --select_alignement SELECT_ALIGNEMENT
                        selection syntax for alignement (backbone)
  -sr SELECT_RMSD, --select_rmsd SELECT_RMSD
                        selection syntax for RMSD
  -m METHOD, --method METHOD
                        method for clustering: single; complete; average;
                        weighted; centroid; median and ward
  -cc CUTOFF, --cutoff CUTOFF
                        cutoff for clusterization from hierarchical clusturing
                        with Scipy. If you choose to clic on the graph, cutoff
                        will be the clicked value otherwise the cutoff will be
                        0.5*max Distance value
  -ng NGROUP, --ngroup NGROUP
                        number of group wanted. Use the maxclust method to
                        clusterize in this case
  -i INTERACTIVE, --interactive INTERACTIVE
                        Interactive mode for distance matrix (Y/n)
```
    
## USAGE : 
```
python TrajectoryClustering.py -f traj.xtc -t TOPOL.pdb -st protein -sr backbone
python TrajectoryClustering.py -f traj.trr -t TOPOL.pdb -sr residue 10 to 30 and chainid 1 -l clustB.log -cc 5 -i n
#For PDB trajectory, d'ont use the **-t** argument
python TrajectoryClustering.py -f traj.pdb -st protein -sr backbone
```

## OUTPUT EXAMPLE
#### Cluster result (structure)
PDBs are saved in a new folder Cluster_PDB.  
You can find in the PDB name the *cluster number*, *size* and the *representative frame*
of the cluster (*ie* the frame number of the saved structure)  
Example: C1-f78-s93.pdb correspond to the cluter 1 made of 93 structures and the
saved frame is the frame 78.

#### Logfile
In the log file you will find all arguments given to the programm with
cluster information:
 - **size**: number structure in the cluster
 - **representative frame**: frame with lowest RMSD between all other frame of the cluster
 - **Members**: all frame belonging to the cluster
 - **spread**: mean RMSD between all frame in the cluster

#### Dendrogram
If you choose to use the matplotlib windows for the final clustering or a cutoff 
(with the **-cc** argument), you will have a picture of the dendrogram.  
The name of this file will be the same as the logfile with a ".png" extension 
example: example.log --> example.png
![alt text](https://github.com/tubiana/TrajectoryClustering/blob/master/examples/example.png "Dendrogram example")  
The grey horizontal line is the cutoff value used

## License
This program is under the GNU GPLv3 license, which means that anyone who 
distributes your code or a derivative work to make the source available under 
the same terms, and also provides an express grant of patent rights from 
contributors to users.

## TODO
 - add quote in commande line writing (in logfile)