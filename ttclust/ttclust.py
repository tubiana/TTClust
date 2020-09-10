#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Thibault TUBIANA"
__license__ = "GNU GPLv3"
__date__ = "2018/02"

import argparse
import datetime
import glob
# ==============================================================================
#                     MODULES
# ==============================================================================
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import operator
import progressbar as pg
import scipy.cluster.hierarchy as sch

try :
    #This if the "builded" import version
    from .version import __version__
except:
    #for usage from sources
    from version import __version__

from prettytable import PrettyTable
from sklearn import manifold

if os.name == "posix":  # Only on linux
    try:
        import argcomplete
    except:
        print("argcomplete not detected, if you want to use autocompletetion")
        print("install argcomplete module.")
        print("See https://github.com/kislyuk/argcomplete for more info")
        pass
# ==============================================================================
#                     GLOBAL VARIABLES
# ==============================================================================
WIDGETS = [pg.Bar('>'), ' ', pg.ETA(), ' ', pg.ReverseBar('<')]
COORDS = []
COLOR_LIST = ["red", "blue", "lime", "yellow",
              "darkorchid", "deepskyblue",
              "orange", "brown", "gray", "black",
              "darkgreen", "navy"]
DPI = 600


# ==============================================================================
#                          CLASS
# ==============================================================================

class Cluster:
    """
    DESCRIPTION
    Simple cluster class object (countains frames numbers, spread, size, ID
    and representative frame)
    """

    def __init__(self, num):
        self.frames = []
        self.spread = -1
        self.size = -1
        self.id = num
        self.representative = -1


# ==============================================================================
#                     TOOL FONCTIONS
# ==============================================================================

def printScreenLogfile(string):
    """
    DESCRIPTION
    Print string on screen and write it on logfile
    Args:
        string (string): string to print and write
    """
    print(string)
    LOGFILE.write("{}\n".format(string))
    LOGFILE.flush()  # forcing the writing by flushing the buffer


def write_command_line():
    """
    DESCRIPTION
    Write commande line with quote for -st, -sr and -sa arguments
    """
    LOGFILE.write("command line       : python ")
    i = 0
    while i < len(sys.argv):
        LOGFILE.write("{} ".format(sys.argv[i]))
        if sys.argv[i] in ["-st", "-sr", "-sa"]:
            i += 1
            LOGFILE.write("\"{}\" ".format(sys.argv[i]))
        i += 1
    LOGFILE.write("\n")


def init_log(args, mdtrajectory):
    """
    DESCRIPTION
    initialyse the logfile with some information
    ----
    Args:
        args (dict): dictionnary of all arguments (argparse)
    """
    topo = args["top"]
    traj = args["traj"]
    selection_string = args["select_traj"]
    select_align = args["select_alignement"]
    select_rmsd = args["select_rmsd"]
    logname = os.path.splitext(args["logfile"])[0]

    LOGFILE.write("========================================================\n")
    LOGFILE.write("====================  TTCLUST {}  ===================\n" \
                  .format(__version__))
    LOGFILE.write("========================================================\n")
    LOGFILE.write("\n")

    LOGFILE.write("************ General information ************\n")
    LOGFILE.write("software version   : {}\n".format(__version__))
    LOGFILE.write("Created on         : {}\n".format(datetime.datetime.now()))
    write_command_line()
    LOGFILE.write("DESTINATION FOLDER : {}\n".format(os.getcwd() + "/" + logname))
    LOGFILE.write("ARGUMENTS : \n")
    LOGFILE.write("  Selection string :\n")
    LOGFILE.write("      Atoms selected in trajectory = {} \n".format(
        selection_string))
    LOGFILE.write("      Atoms selected for alignement = {} \n".format(
        select_align))
    LOGFILE.write("      Atoms selected for RMSD = {} \n".format(select_rmsd))
    LOGFILE.write("  trajectory file  : {} \n".format(','.join(traj)))
    LOGFILE.write("   Stride          : {} \n".format(args["stride"]))
    LOGFILE.write("   Number of frames  : {} \n".format(mdtrajectory.n_frames))
    LOGFILE.write("   Number of atoms  : {} \n".format(mdtrajectory.n_atoms))
    LOGFILE.write("  topology file    : {} \n".format(topo))
    LOGFILE.write("  method used of clusterring : {}".format(args["method"]))
    LOGFILE.write("\n\n")
    if args["ngroup"]:
        LOGFILE.write("  Number of cluster asked: {}\n".format(args["ngroup"]))
    if args["cutoff"]:
        LOGFILE.write("  cutoff for dendrogram clustering: {}\n".format("cutoff"))


def extract_selected_atoms(selection, traj, logname, save=False):
    """
    DESCRIPTION
    Return a trajectory with only atoms given in arguments (trough the
    selection string)
    eg: if you only want the trajectory of the chain A you can call
    traj_A = extract_selected_atoms("chainid 0", traj)
    ----
    Args:
        selection (string): selection string (mdtraj syntax)
        traj (mdtraj.trajectory): initial trajectory
        logname (string): logname (basename) for folder location.
    Returns:
        subtraj (mdtraj.trajecotry): subtrajectory of selected atoms
    """
    try:
        subtraj = traj.atom_slice(traj.top.select(selection))
        subtraj.center_coordinates()
        if save:
            printScreenLogfile("NOTE : 'st' argument given. I will save the subtrajectory"
                               " in {0}/{0}.xtc and topology file as {0}/{0}.pdb".format(logname))
            subtraj[0].save_pdb("{0}/{0}.pdb".format(logname))
            subtraj.save_xtc("{0}/{0}.xtc".format(logname))
        return subtraj
    except:
        print("ERROR : there is an error with your selection string")
        print("        SELECTION STRING : ")
        print(("        {}".format(selection)))
        print("        > Please check 'http://mdtraj.o"
              "rg/latest/atom_selection.html'")
        exit(1)


def send_error_message(calc_type, selection_string, other=""):
    """
    DESCRIPTION
    Print information regarding your selection string if this one is not
    recognized by mdtraj
    Args:
        calc_type
    """
    print(("ERROR : {} selection string not valid".format(calc_type)))
    print(("        >{}".format(selection_string)))
    if not other == "":
        print("        >{}".format(other))
    exit(1)


def improve_nucleic_acid(selection_string):
    """
    DESCRIPTION
    improve RNA and DNA selection
    ---
    Args:
        selection_string (string) : selection string (for MDtraj)
    Return:
        selection_string (string) : improved RNA selection string
    """
    dna = "(resname =~ '(5|3)?D([ATGC]){1}(3|5)?$')"
    rna = "(resname =~ '(3|5)?R?([AUGC]){1}(3|5)?$')"
    backbone_na = "rna or dna and (name =~ \"(P)|(O[35]')|(C[3-5]')\")"
    base = "(rna or dna) and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H"
    base_rna = "(rna) and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H"
    base_dna = "dna and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H"
    if 'base_rna' in selection_string:
        selection_string = selection_string.replace('base_rna', base_rna)
    if 'base_dna' in selection_string:
        selection_string = selection_string.replace('base_dna', base_dna)
    if 'base' in selection_string:
        selection_string = selection_string.replace('base', base)
    if 'backbone_na' in selection_string:
        selection_string = selection_string.replace('backbone_na', backbone_na)
    if 'dna' in selection_string:
        selection_string = selection_string.replace('dna', dna)
    if 'rna' in selection_string:
        selection_string = selection_string.replace('rna', rna)
    return selection_string


def return_selection_atom(use_for, traj, args):
    """
    DESCRIPTION
    return indices of selected atoms.
    ----
    Args:
        use_for (string): witch selection string was wrong ? (sr/sa/sr)
        traj (mdtraj.trajectory): trajectory
        args (dictionnary): dictionnary of all arguments
    return:
        selection (array): array of selected atoms (atom numbers in index 0 based).
    """
    selection_string = args["select_alignement"]
    try:
        selection = traj.top.select(selection_string)
    except ValueError:
        send_error_message(use_for, selection_string, "Keyword not recognize")

    if len(selection) == 0:
        if selection_string == "backbone":
            selection = traj.top.select(improve_nucleic_acid("backbone_na"))
            args["select_alignement"] = "backbone_na"
            args["select_rmsd"] = "backbone_na"
        # If DNA or RNA wasn't selected, stop the program.
        if len(selection) == 0 or selection == None:
            send_error_message(use_for, selection_string, "Selection list EMPTY")
        else:
            print("NOTE : Nucleic acids found.")
            print("       Automatic switch to nucleic acid mode")

    return selection


def save_dist_mat(distmat, rmsd_string):
    """
    DESCRIPTION
    Save the numpy matrix to reused afterward
    ----
    Args:
        distmat (numpy matrix): distance matrix
        rmsd_string (str) : selection string for rmsd calculation
    return:
        None
    """
    if rmsd_string:
        name = rmsd_string.replace(" ", "_")
    else:
        name = "matrix_all"
    np.save(name, distmat)
    printScreenLogfile("Saving distance matrix : {0}.npy".format(name))


def reorder_cluster(clusters):
    """
    DESCRIPTION
    Reorder the clusters number to have the first frame belonging to the
    cluster 1.
    ---
    Args :
        Clusters_labels(list): list of clusters label.
    """
    dict_order = {}
    for i in range(len(clusters)):
        dict_order[clusters[i].id] = clusters[i].frames[0]
    # Evaluate order
    sorted_clusters = sorted(dict_order.items(), key=operator.itemgetter(1))
    for i in range(len(sorted_clusters)):
        dict_order[sorted_clusters[i][0]] = i + 1  # i+1 == reorder cluster number
    # reordering
    for i in range(len(clusters)):
        clusters[i].id = dict_order[clusters[i].id]


# ==============================================================================
#                     FONCTIONS
# ==============================================================================

# @Gooey(progress_regex=r"\|( |>)*\| ETA: +([0-9]|-)+:([0-9]|-)+:([0-9]|-)+ \|( |<)*\|",
# disable_progress_bar_animation=True)
def parseArg():
    """
    This fonction will the list of pdb files and the distance
    @return: dictionnary of arguments
    Ex :
    python Cluster_Analysis.py -f *.pdb -s A:1-30:CA
    """
    arguments = argparse.ArgumentParser(description="This program was developped in order to clusterize "
                                                    "molecular dynamictrajectories (Amber, gromacs, chamm, namd, PDB)")
    try:
        argcomplete.autocomplete(arguments)
    except:
        pass
    arguments.add_argument('-f', "--traj", help="trajectory file(s). You can give a list of files.", required=True, nargs='+')
    arguments.add_argument('-t', '--top', help="topfile", default=None)
    arguments.add_argument('-s', '--stride', help="stride (read every Xth frames", default=1,type = int)
    arguments.add_argument('-l', '--logfile', help="logfile (default : clustering.log). The "
                                                   "name of your output file will be the basename (name before the extention "
                                                   "of this logfile", default="clustering")
    arguments.add_argument('-st', '--select_traj', help="selection syntax for "
                                                        "Don't forget to add QUOTES besite this selection string."
                                                        "trajectory extraction (default : all).", default="all")
    arguments.add_argument('-sa', '--select_alignement', help="selection syntax"
                                                              " for alignement (default : backbone). Don't forget to add QUOTES besite this "
                                                              "selection string."
                                                              " If you don't want aligment use \"none\".",
                           default="backbone")
    arguments.add_argument('-sr', '--select_rmsd', help="selection syntax for "
                                                        " RMSD (default : backbone). Don't forget to add QUOTES "
                                                        "besite this selection string.", default="backbone")

    # Clustering arguments
    arguments.add_argument('-m', '--method', help="method for clustering : single "
                                                  "; complete; average; weighted; centroid; median. (ward)",
                           default="ward")
    arguments.add_argument('-cc', "--cutoff", help="cutoff for clusterization from "
                                                   "hierarchical clusturing with Scipy", default=None)
    arguments.add_argument('-ng', "--ngroup", help="number of group asked. Use the "
                                                   "maxclust method to clusterize in this case", default=None)
    arguments.add_argument('-aa', "--autoclust", help="Autoclustering (Y/n)", default="Y")

    # Interactive mode for distance matrix:
    arguments.add_argument('-i', '--interactive', help="Use saved distance matrix ? (Y/n)", default="Y")

    args = vars(arguments.parse_args())

    # Activate autoclustering if autoclust is True and if no values was specified for ngroup or cutoff
    if args["autoclust"] in ["Y", "y"]:
        args["autoclust"] = True
    else:
        args["autoclust"] = False

    if (args["autoclust"] == True) and (args["ngroup"] == None) and (args["cutoff"] == None):
        args["ngroup"] = "auto"

    return (args)


def ask_choice(args, name):
    """
    DESCRIPTION
    If a distance matrix file is found (the name of the matrixe is the same
    as the rmsd selection string), the programe ask the the user if he want to
    use it
    ---
    Args:
        args (dict): all arguments in a dictionary
        name (string): file name
    Return:
        name (string): if we use the matrix file we send back the file name
        None (None): otherwise we send back nothing
    """
    if args["interactive"].upper() == "Y":
        print("Interactive mode disabled. I will use the saved matrix")
        return name

    print(" I found a distance matrix ({0}) saved. Do you want to use it ?".format(name))
    print("    y/Y - YES")
    print("    n/N - NO")
    print("    o/O - find all other .npy distance matrix")
    # aks the user what to do!
    if sys.version_info[0] == 2:
        choice = raw_input()
    else:
        choice = input()
    # evaluate answer.
    if choice.upper() == "Y":  # I want to use it!
        printScreenLogfile(" >Distance matrix file detected : {0}".format(name))
        return (name)
    elif choice.upper() == "N":  # don't want to use it.. Recalculate!
        print("Calculation mode activated")
        return None
    elif choice.upper() == "O":  # I want to use another npy distance matrix
        npy_files = glob.glob("*.npy")
        for i, file in enumerate(npy_files):
            print("  {0} - {1}".format(i + 1, file))
        print(" -->Please chooce and press Enter")
        # Check if the user give a good answer
        choice_file = input()
        try:
            name = npy_files[int(choice_file) - 1]
            return name
        except:
            print("I didn't understand. Please try again")
            print("........")
            return ask_choice(args, name)
    else:
        print("I didn't understand. Please try again")
        print("........")
        return ask_choice(args, name)


def search_dist_mat(rmsd_string, args):
    """
    Search if the distance matrix already exist
    ----
    Args:
        rmsd_string (str) : name of the numpy matrix
    """
    if rmsd_string:
        name = rmsd_string.replace(" ", "_")
    else:
        name = "matrix_all"
    # Searching all npy file in the folder
    npy_files = glob.glob("*.npy")
    if not name[-4:] == ".npy":
        name += ".npy"

    if name in npy_files and not args["interactive"].lower() == "n":
        return ask_choice(args, name)
    else:
        return None


def calculate_representative_frame_spread(clusters_list, DM):
    """
    DESCRIPTION
    Choose the representative frame by calculating the mean RMSD of each
    structures of the cluster agains the others
    ----
    Args:
        clusters_list (array): list of all clusters
        DM (Numpy matrix): distance matrix for each frames
    """
    print("Searching for representative frames")

    for n, cluster in enumerate(clusters_list):
        frames = cluster.frames
        mean_rmsd_per_frame = {}
        # first loop  : first frame
        for frame_i in frames:
            mean_rmsd_per_frame[frame_i] = 0
            # we will add the rmsd between theses 2 frames and then calcul the
            # mean
            for frame_j in frames:
                # We don't want to calcul the same frame.
                if not frame_j == frame_i:
                    # we add to the corresponding value in the list of all rmsd
                    # the RMSD betwween frame_i and frame_j
                    mean_rmsd_per_frame[frame_i] += DM[frame_i - 1, frame_j - 1]
            # mean calculation
            mean_rmsd_per_frame[frame_i] /= (len(frames))

            # Representative frame = frame with lower RMSD between all other
            # frame of the cluster
            repre = min(mean_rmsd_per_frame, key=mean_rmsd_per_frame.get)
            cluster.representative = repre

            # spread = mean rmsd in all the cluster (*10 to have angstöm)
            cluster.spread = sum(mean_rmsd_per_frame.values()) / len(frames)
            cluster.spread *= 10


def create_DM(traj, args):
    """
    DESCRIPTION
    Calcul the distance matrix
    ---
    Args:
        traj (mdtraj.trajectory): trajectory
        args (dict): all arguments in dictionary
    return:
        distances (numpy matrix): distance matrix
    """
    # Get Atoms indices from selection string

    if args["select_alignement"] != "none":
        alignement_selection = return_selection_atom(use_for="ALIGNEMENT",
                                                     traj=traj,
                                                     args=args)

        # Trajectory superposition  (aligment)
        traj_aligned = traj.superpose(traj[0],
                                      atom_indices=alignement_selection,
                                      parallel=True)
    else:
        traj_aligned = traj

    select_align = improve_nucleic_acid(args["select_alignement"])
    untouch_rmsd_string = args["select_rmsd"]
    rmsd_string = improve_nucleic_acid(args["select_rmsd"])

    if rmsd_string:
        print("NOTE : Extraction of subtrajectory for time optimisation")
        traj_aligned = extract_selected_atoms(rmsd_string, traj_aligned, args["logname"])
    # matrix initialization
    distances = np.empty((traj.n_frames, traj.n_frames))

    # Searching if a distance file already exist
    distance_file = search_dist_mat(untouch_rmsd_string, args)

    # If a distance matrix file was found and choosed, we load it.
    if distance_file:
        printScreenLogfile(" >Distance Matrix File Loaded!")
        return np.load(distance_file)
    else:  # otherwise
        pbar = pg.ProgressBar(widgets=WIDGETS, maxval=traj.n_frames).start()
        counter = 0
        # Pairwise RMSD calculation (matrix n²)
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj_aligned, traj_aligned, frame=i)
            pbar.update(counter)
            counter += 1
        pbar.finish()

        # Finaly, we save the matrix if we want to load it again afterward
        print("Calculation ended - saving distance matrix")
        save_dist_mat(distances, args["select_rmsd"])

        return distances


def onclick(event):
    """
    DESCRIPTION
    This function is used to get coordinate of the mouse on the matplotlib
    windows.
    """
    ix, iy = event.xdata, event.ydata

    global COORDS
    COORDS.append((ix, iy))

    # juste one clic
    if len(COORDS) == 1:
        plt.close(1)


def return_mapping_cluster(labels):
    """
    DESCRIPTION
    assign cluster to frames.
    ---
    Args:
        labels (list) : each "label" is the cluster number
    Returns:
        clusters_list : list of clusters
    """
    # cluster_with_frame_number = defaultdict(lambda : [])
    clusters_list = []
    for cluster_num in set(labels):
        clusters_list.append(Cluster(cluster_num))  # create new instance of cluster

    for i, cluster_num in enumerate(labels):
        clusters_list[cluster_num - 1].frames.append(i)
        # for DEBUG
        if cluster_num != clusters_list[cluster_num - 1].id:
            print ("{0} - {0}".format(cluster_num, clusters_list[cluster_num - 1]))
            sys.exit(1)

    for cluster in clusters_list:
        cluster.size = len(cluster.frames)
    # print(mapping)
    return clusters_list


def segments_gain(p1, v, p2):
    # From https://datascience.stackexchange.com/questions/6508/k-means-incoherent-behaviour-choosing-k-with-elbow-method-bic-variance-explain
    vp1 = np.linalg.norm(p1 - v)
    vp2 = np.linalg.norm(p2 - v)
    p1p2 = np.linalg.norm(p1 - p2)
    return np.arccos((vp1 ** 2 + vp2 ** 2 - p1p2 ** 2) / (2 * vp1 * vp2)) / np.pi


def auto_clustering(matrix):
    """
    DESCRIPTION
    Autoclustering function based on sklearn (for now) and the elbow method
    Based on: https://datascience.stackexchange.com/questions/6508/k-means-incoherent-behaviour-choosing-k-with-elbow-method-bic-variance-explain
    """
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist

    distorsions = []
    K = range(2, 15)
    for k in K:
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(matrix)

        distorsions.append(sum(np.min(cdist(matrix, kmeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0])

    criterion = np.array(distorsions)
    criterion = (criterion - criterion.min()) / (criterion.max() - criterion.min())

    # Compute the angles
    seg_gains = np.array([0, ] + [segments_gain(*
                                                [np.array([K[j], criterion[j]]) for j in range(i - 1, i + 2)]
                                                ) for i in range(len(K) - 2)] + [np.nan, ])

    # Get the first index satisfying the threshold
    seg_threshold = 0.99  # Set this to your desired target

    kIdx = np.argmax(seg_gains > seg_threshold)

    kmeans = KMeans(n_clusters=kIdx)
    kmeans.fit(matrix)
    # return(labels)
    return (kIdx)


def create_cluster_table(traj, args):
    """
    DESCRIPTION
    Clustering function!
    Create a list with the cluster number of earch frame
    eg: [1,2,1,2,2,2,1] which mean that:
        the cluster 1 is composed of the frames 1,3,7
        the cluster 2 is composed of the frames 2,4,5,6
    Args:
        traj (mdtraj.trajectorie): trajectory file
        args (dict): all arguments in a dictionary
    Return:
        Distances (numpy matrix): Distance matrix
        clustering_result (list): cluster unmber list for each frame (index)
    """
    print("         creating distance matrix")
    distances = create_DM(traj, args)

    select_align = args["select_alignement"]
    untouch_select_rmsd = args["select_rmsd"]
    select_rmsd = improve_nucleic_acid(args["select_rmsd"])
    cutoff = args["cutoff"]
    ncluster = args["ngroup"]
    # Creation of the distance matrix

    if select_rmsd == None:
        select_rmsd = "None"

    linkage_file = search_dist_mat(untouch_select_rmsd + " linkage " + args["method"], args)

    if linkage_file:
        linkage = np.load(linkage_file)
    else:
        print("         Matrix shape: {}".format(distances.shape))
        print("         Scipy linkage in progress. Please wait. It can be long")
        try:
            linkage = sch.linkage(distances, method=args["method"])
        except:
            printScreenLogfile("ERROR : method name given for clustering didn't recognized")
            printScreenLogfile("      : methods are : single; complete; average; weighted; centroid; ward.")
            printScreenLogfile("      : check https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/sc"
                               "ipy.cluster.hierarchy.linkage.html for more info")
            sys.exit(1)
        print("         >Done!")
        print("         ...Saving linkage matrix...")
        save_dist_mat(linkage, untouch_select_rmsd + " linkage " + args["method"])
        print("         >Done!")

    # if a cuttof for distance cuting is given
    if cutoff:
        clustering_result = sch.fcluster(linkage, cutoff, "distance")
    # otherwise we choose it on the screen by cliking on the matplotlib windows
    # If a number of wanted cluster is given
    elif ncluster:
        if ncluster == "auto":
            # clustering_result = sch.fcluster(linkage,4, 'distance')
            # clustering_result = auto_clustering(distances)
            ncluster = auto_clustering(distances)
        clustering_result = sch.fcluster(linkage, t=ncluster, criterion="maxclust")
        # print(len(clustering_result))
        n_group = len(np.unique(clustering_result))
        cutoff = linkage[-(n_group - 1), 2]
    else:
        clicked = False
        while not clicked:
            fig = plt.figure()
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.title("Please click where you want to build clusters")
            sch.dendrogram(linkage)
            plt.show()
            try:
                cutoff = COORDS[0][1]
                clicked = True
                clustering_result = sch.fcluster(linkage, cutoff, "distance")
            except:
                print("ERROR : PLEASE CLICK ON THE DENDROGRAM TO CHOOSE YOUR CUTOFF VALUE")
                print("        Please wait during the dendrogram regeneration")
            plt.close()
    printScreenLogfile("  cutoff for clustering : {:.2f}".format(float(cutoff)))
    return distances, clustering_result, linkage, cutoff


def write_representative_frame(traj, cluster, logname):
    """
    DESCRIPTION
    Write representative frame of a cluster
    ----
    Args:
        traj (mdtraj.trajectory): trajectory
        cluster (Cluster): a Cluster object
        logname (string): output logfile
    """
    cluster_num = cluster.id
    frame = cluster.representative
    size = cluster.size
    #bugfix in 4.6.8
    traj[frame].save_pdb("{}/C{}-f{}-s{}.pdb".format(logname,
                                                     cluster_num,
                                                     frame+1, # +1 to get the 1 index based frame
                                                     size))


def get_cmap(num_cluster):
    """
    DESCRIPTION
    This function will return cmap that will be used for graphics
    Args:
        num_cluster (int) : number of cluster
    Returns:
        cmap (matplotlib cmap) : matplotlib color map
    """
    global COLOR_LIST
    # DEFINE COLOR MAP
    # if two much cluster number to define colours by hand

    if num_cluster > len(COLOR_LIST):
        print()
        cmap = "rainbow_r"
    else:
        # imshow take the last color for the last group (if 3 cluster, color of
        # clusters 3 will be brown")
        COLOR_LIST = COLOR_LIST[:num_cluster]
        cmap = mpl.colors.ListedColormap(COLOR_LIST)
    return cmap


def plot_barplot(clusters_list, logname, size, traj):
    """
    DESCRIPTION
    This function is used to plot the linear barplots.
    Args:
        clusters_list (list) : list of cluster label in order or appearance
        logname (str) : output logname
        size (int): number of frames
    Returns:
        colors_list (list of list) : list of colors in RGBA format
    """
    # order clusters_labels by order of appearance in the trajectory
    clusters_number_ordered = [0] * size
    # Order clusters_labels by cluster order.
    for cluster in clusters_list:
        for frame in cluster.frames:
            clusters_number_ordered[frame] = cluster.id

    # DEFINE COLOR MAP
    cmap = get_cmap(len(clusters_list))

    data = np.asmatrix(clusters_number_ordered)

    fig, ax = plt.subplots(figsize=(10, 1.5))

    # get time if present
    if traj.time.sum() < 0.0000005:
        timeMin, timeMax = 0, np.shape(data)[1]
        plt.xlabel("Frame")

    else:
        try:
            timeMin, timeMax = traj.time[0] / 1000, traj.time[-1] / 1000  # For time in ns.
            plt.xlabel("Time (ns)")
            if timeMax >= 1000:
                timeMin = timeMin / 1000
                timeMax = timeMax / 1000
                plt.xlabel("Time ($\mu$s)")
        except:
            timeMin, timeMax = 0, np.shape(data)[1]

    im = plt.imshow(data, aspect='auto', interpolation='none', cmap=cmap, extent=[timeMin, timeMax, 1, 0])

    plt.tight_layout()
    plt.tick_params(axis="y",
                    which='both',
                    left=False,
                    right=False,
                    labelleft=False)
    plt.tick_params(axis="x",
                    direction="out",
                    which='both',
                    top=False)
    colors_list = (im.cmap(im.norm(np.unique(clusters_number_ordered))))

    plt.savefig("{0}/{1}-linear.png".format(logname,
                                            logname.split(os.sep)[-1]),
                dpi=DPI,
                transparent=True)
    plt.close()
    return colors_list


def plot_hist(clusters_list, logname, colors_list):
    """
    DESCRIPTION
    This function is used to plot a histogram with the cluster size.
    Args:
        clusters_list (list): list of cluster label in order or appearance
        logname (str): output logname
        colors_list (list): list of colors
    Returns:
        None
    """
    if mpl.__version__[0] == "2":
        STYLE = "classic"
        if STYLE in plt.style.available:
            plt.style.use(STYLE)
    values = []
    labels = []
    for cl in clusters_list:
        # occurence.append((cl.id, cl.size))
        values.append(cl.size)
        labels.append(cl.id)
    # Sorting occurence dict by cluster size

    #### Configuration plot
    width = 0.7  # bars size
    index = np.arange(len(values))  # the x locations for the groups
    fig, ax = plt.subplots()

    bp = ax.bar(index, values, width, color=colors_list, label="Cluster size")
    # add value on top of bars, adapted from matplotlib doc
    for rect in bp:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.0 * height,
                '%d' % int(height),
                ha='center', va='bottom')

    plt.xlabel("Clusters")
    plt.ylabel("Number of members")
    plt.title("Distribution within clusters")
    plt.xticks(index + (width / 2), labels)
    plt.tight_layout()

    plt.savefig("{0}/{1}-hist.png".format(logname,
                                          logname.split(os.sep)[-1]),
                dpi=DPI, transparent=True)
    plt.close()


def plot_distmat(distances, logname):
    """
    DESCRIPTION
    
    """

    fig, ax = plt.subplots()

    plt.imshow(distances, interpolation='none', origin='lower')
    plt.colorbar()
    plt.xlabel("Frame")
    plt.ylabel("Frame")
    plt.title("RMSD between frames (nm)")
    plt.savefig("{0}/{1}-distmat.png".format(logname,
                                             logname.split(os.sep)[-1]),
                dpi=DPI,
                transparent=True)
    plt.close()


def plot_dendro(linkage, logname, cutoff, color_list, clusters_list):
    """
    DESCRIPTION
    This function will create the dendrogram graph with the corresponding
    cluster color.
    Args:
        linkage (numpy array) : linkage matrix
        logname (str) : output logfile name
        cutoff (float) : cutoff used for clustering
        color_list (list) : HEX code color for each cluster
        clusters_list (list) : list of all cluster (Cluster object)
    Returns:
        None
    """
    if mpl.__version__[0] == "2":
        STYLE = "classic"
        if STYLE in plt.style.available:
            plt.style.use(STYLE)
    fig = plt.figure()
    # Convert RGB color to HEX color
    color_hex = [mpl.colors.rgb2hex(x) for x in color_list]
    sch.set_link_color_palette(color_hex)
    # clusters_list
    color_member = {}
    for cl in clusters_list:
        for frm in cl.frames:
            color_member[frm] = mpl.colors.rgb2hex(color_list[cl.id - 1])

    # Attribute the correct color for each branch.
    # adapte from Ulrich Stern code in StackOverflow http://stackoverflow.com/a/38208611
    link_cols = {}
    for i, i12 in enumerate(linkage[:, :2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(linkage) else color_member[x] for x in i12)
        link_cols[i + 1 + len(linkage)] = c1 if c1 == c2 else "#808080"

    # Dendrogram creation
    # Override the default linewidth.
    den = sch.dendrogram(linkage, color_threshold=float(cutoff), above_threshold_color="#808080",
                         link_color_func=lambda x: link_cols[x])

    # Graph parameters
    plt.title("Clustering Dendrogram")
    ax = plt.axes()
    ax.set_xticklabels([])
    plt.axhline(y=float(cutoff), color="grey")  # cutoff value vertical line
    ax.set_ylabel("Distance (AU)")
    ax.set_xlabel("Frames")

    plt.savefig("{0}/{1}-den.png".format(logname,
                                         logname.split(os.sep)[-1]),
                format="png", dpi=DPI,
                transparent=True)
    plt.close()


def symmetrize_matrix(matrix):
    """
    DESCRIPTION
    This function will make a symmetric matrix from another matrix (from
    the top part of the matrix). Usefull of multidimentional scaling.
    Args:
        matrix (np.array)
    Returns:
        matrix_sym (np.array)
    """
    dim = matrix.shape[0]
    matrix_sym = np.copy(matrix)
    for i in range(dim):
        for j in range(i, dim):
            matrix_sym[j, i] = matrix[i, j]
    return matrix_sym


def plot_2D_distance_projection(rmsd_m, clusters_list, colors, logname):
    """
    DESCRIPTION
    This function will create a 2D distance projection graph with the MDS methods
    Args:
        rmsd_m (np.array) : rmsd matrix (between clusters)
        clusters_list (list of Cluster): list of Clusters
    Return:
        None
    """
    labels = range(1, len(clusters_list) + 1)
    # 1 - value normalisation (make value between 0 and 1) of RMSD matrix
    rmsd_norm = rmsd_m / np.max(rmsd_m)
    symmetrize_matrix(rmsd_norm)

    rmsd_norm = symmetrize_matrix(rmsd_norm)
    # 2 - create the MDS methods
    # mds = manifold.MDS(n_components=2, dissimilarity="euclidean", random_state=4)
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed")  # , random_state=2)

    # 3 - MDS projection
    rmsd_mds = mds.fit(rmsd_norm)
    # rmsd_mds = mds.fit(rmsd_m)

    # 4 - get X/Y coords
    coords = rmsd_mds.embedding_

    # 5 - get spread and normalyse
    spreads = []
    for clust in clusters_list:
        spreads.append(clust.spread)
    spreads = np.array(spreads)
    #    spreads_norm = spreads / np.max(spreads)
    # minspread = np.min(spreads_norm)+0.05*np.min(spreads_norm)
    radii = np.pi * (25 * (spreads) ** 2)  # radii = 5 to 20
    x = coords[:, 0]
    y = coords[:, 1]

    # 6 - plot graph
    fig = plt.figure()
    ax = plt.subplot(111)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    scatter = ax.scatter(x, y, s=radii, c=colors, alpha=0.5)
    for label, x, y in zip(labels, x, y):
        plt.annotate(label,
                     xy=(x, y),
                     ha='left', va='bottom', fontsize=8)

    # set the same axis for X and Y

    lims = []
    lims.extend(ax.get_xlim())
    lims.extend(ax.get_ylim())
    ax.set_ylim((min(lims), max(lims)))
    ax.set_xlim((min(lims), max(lims)))

    plt.title("Relative distance between clusters")
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labelbottom='off')  # labels along the bottom edge are off
    plt.tick_params(
        axis='y',  # changes apply to the y-axis
        which='both',  # both major and minor ticks are affected
        left="off",
        right="off",
        labelleft='off')  # labels along the bottom edge are off

    # 7 - circle bar
    max_size = max(radii)
    min_size = min(radii)
    min_color = colors[np.argmin(radii)]
    max_color = colors[np.argmax(radii)]

    # add transparency
    min_color[-1] = 0.5
    max_color[-1] = 0.5
    leg_min = plt.scatter([], [], s=min_size, edgecolor='black', color=min_color)
    leg_max = plt.scatter([], [], s=max_size, edgecolor='black', color=max_color)
    labels = ["{:.2f}".format(min(spreads)), "{:.2f}".format(max(spreads))]

    legend = ax.legend([leg_min, leg_max], labels,
                       ncol=1, frameon=False, fontsize=8,
                       handlelength=2, loc="upper right", borderpad=1.8, handletextpad=1,
                       scatterpoints=1, bbox_to_anchor=(1.3, 0.9))
    legend.set_title('Spread radius', prop={"size": "small"})

    # Add Text for distance information
    min_rmsd = np.min(rmsd_m[np.nonzero(rmsd_m)])
    max_rmsd = np.max(rmsd_m[np.nonzero(rmsd_m)])
    text_distance = ("RMSD\n   min : {:.2f}$ \AA$\n   max : {:.2f} $\AA$".format(min_rmsd, max_rmsd))

    # plt.gca().add_artist(legend1)
    ax.annotate(text_distance, xy=(1.05, 0.5), xycoords="axes fraction", fontsize="small")

    plt.savefig("{0}/{1}-dist.png".format(logname,
                                          logname.split(os.sep)[-1]),
                format="png",
                dpi=DPI, transparent=True)
    plt.close()


def generate_graphs(clusters_list, output, size, linkage, cutoff, distances, traj):
    """
    DESCRIPTION
    Create a linear cluster mapping graph where every frame is printed as a
    colored barplot
    Args:
        clusters_list (list): list of cluster
        output (string): output name for graph
        size (int): number of frames
        linkage (numpy array): matrix linkage
        cutoff (float): cutoff distance value for clustering (in the dendogram)
        distances(numpy array): distance matrix
        traj (Trajectory): trajectory for time usage in axis barplot
    Return:
        colors_list (list) to be used with 2D distance projection graph
    """
    colors_list = plot_barplot(clusters_list, output, size, traj)
    plot_dendro(linkage, output, cutoff, colors_list, clusters_list)
    plot_hist(clusters_list, output, colors_list)
    if (distances.shape[0] < 10002):
        plot_distmat(distances, output)
    else:
        printScreenLogfile("Too many frames (>=10002)! The RMSD distance matrix will not be generated")
    return colors_list


def get_RMSD_cross_cluster(clusters_list, distances, logname):
    """
    DESCRIPTION
    This function will get the RMSD between all representativ frames of all
    clusters. Print it to the console and write it on the logfile
    Args:
        clusters_list (list): list of all clusters
        distances (np matrix): rmsd matrix
        logname (string): logname (same as logfile)
    returns:
        RMSD_matrix (np.array) : RMSD matrix between clusters
    """
    # 1 - table preparation
    table = PrettyTable()
    n_clusters = len(clusters_list)
    field_names = ["Clusters"] + ["C" + str(x) for x in range(1, n_clusters + 1)]
    table.field_names = field_names
    table.float_format = ".2"
    # + variable which countains all value except 0.0 values (for
    # average calculation)
    non_diag_value = []

    # 2 - RMSD "calculation"
    RMSD_matrix = np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        repr1 = clusters_list[i].representative
        for j in range(i + 1, n_clusters):
            repr2 = clusters_list[j].representative
            # reprx-1 to correspond with the numpy index
            rmsd = distances[repr1][repr2] * 10
            RMSD_matrix[i][j] = RMSD_matrix[j][i] = rmsd
            non_diag_value.append(rmsd)

    # 3 - PrettyTable creation
    for i in range(n_clusters):
        table.add_row(["C" + str(i + 1)] + RMSD_matrix[i].tolist())

    # 4 - print table
    printScreenLogfile("----------------------------")
    printScreenLogfile("RMSD MATRIX BETWEEN CLUSTERS")
    printScreenLogfile(table)
    printScreenLogfile("\nAVERAGE RSMD BETWEEN CLUSTERS : {:.2f}".format(
        np.mean(non_diag_value)))
    np.savetxt("{0}/RMSD_between_clusters.csv".format(logname), RMSD_matrix, delimiter=";")
    return RMSD_matrix


def Cluster_analysis_call(args):
    """
    DESCRIPTION
    Main function of the program : call other function as a pipeline
    Args:
        args (dict): all arguments in a dictionary
    Return:
        traj (Trajectory): simulation trajectory
    """
    trajfile = args["traj"]
    topfile = args["top"]
    select_traj = improve_nucleic_acid(args["select_traj"])
    # Check if "logfile finish with ".log"

    logname = os.path.splitext(args["logfile"])[0]

    # IMPROVE DNA AND RNA SELECTION

    args["select_traj"] = improve_nucleic_acid(args["select_traj"])
    args["select_alignement"] = improve_nucleic_acid(args["select_alignement"])

    print("======= TRAJECTORY READING =======")
    if len(trajfile) == 1:
        trajfile = trajfile[0]
        if topfile == None and trajfile[-4:] == ".pdb":
            traj = md.load_pdb(trajfile)
        else:

            traj = md.load(trajfile,
                           top=topfile, stride=args["stride"])
    elif len(trajfile) > 1:
        print(">Several trajectories given. Will concatenate them.")
        trajList = []
        for t in trajfile:
            if topfile == None and t[-4:] == ".pdb":
                trajList.append(md.load_pdb(t))
            else:
                trajList.append(md.load(t,
                               top=topfile,stride=args["stride"]))
        traj = md.join(trajList)

        #resting the timetable
        if traj.timestep > 0:
            traj.time = np.asarray(list(range(0, int(len(traj) * traj.timestep), int(traj.timestep))))

    else:
        print("ERROR: no trajectory given. TTClust will stop")
        sys.exit(1)

    init_log(args, traj)

    if not select_traj == "all":
        print("======= EXTRACTION OF SELECTED ATOMS =======")
        traj = extract_selected_atoms(select_traj, traj, args["logname"], save=True)

    print("====== Clustering ========")
    distances, clusters_labels, linkage, cutoff = create_cluster_table(traj, args)

    printScreenLogfile("\n**** Cluster Results")
    clusters_list = return_mapping_cluster(clusters_labels)

    print("====== Reordering clusters ======")
    reorder_cluster(clusters_list)
    # reordering the list by the cluster number
    clusters_list.sort(key=operator.attrgetter("id"))
    print("====== Generating Graph ======")
    colors_list = generate_graphs(clusters_list, logname, traj.n_frames, linkage, cutoff, distances, traj)
    print("====== Calc. repr. frame  ======")
    calculate_representative_frame_spread(clusters_list, distances)

    for cluster in clusters_list:
        printScreenLogfile("cluster {}".format(cluster.id))
        printScreenLogfile("    size = {}".format(cluster.size))
        printScreenLogfile("    representative frame={}".format(
            cluster.representative))
        printScreenLogfile("    spread  : {0:.2f} ".format(cluster.spread))
        printScreenLogfile("    Frames : {} ".format(str([x + 1 for x in cluster.frames])))
        write_representative_frame(traj, cluster, logname)

    RMSD_matrix = get_RMSD_cross_cluster(clusters_list, distances, logname)

    plot_2D_distance_projection(RMSD_matrix, clusters_list, colors_list, logname)
    # return trajectory for usage afterwards.
    return traj


def define_LOGFILE(log):
    """
    Define LOGFILE if called from GUI
    """
    global LOGFILE
    LOGFILE = log


def main():
    """
    Execute TTclust

    Return:
        traj (Trajectory): simulation trajectory
    """

    print("********************************************************")
    print("******************  TTCLUST {} *********************".format(
        __version__))
    print("********************************************************")
    print("")

    # We get all arguments
    args = parseArg()
    global LOGFILE

    # add ".log" if the logfile doesn't have extension
    if os.path.splitext(args["logfile"])[1] == "":
        args["logfile"] = args["logfile"] + ".log"

    # create a folder based on the logfile name and write everything inside
    logname = os.path.splitext(args["logfile"])[0]
    args["logname"] = logname
    if not os.path.exists(logname):
        os.makedirs(logname)
    elif not os.path.isdir(logname):  # If a file exist with the same foldername
        os.rename(logname, logname + ".bak")
        print("NOTE : A file with the same folder name was found and rename "
              "into {}.bak".format(logname))
        os.makedirs(logname)

    filename = args["logfile"].split(os.sep)[-1]
    LOGFILE = open("{0}/{1}".format(logname, filename), "w")
    traj = Cluster_analysis_call(args)
    LOGFILE.close()

    # return traj for usage afterwards
    return traj


###############################################################################
#####                               MAIN                                 ######
###############################################################################
if __name__ == "__main__":
    main()  # keep trajectory for usage afterwards (in shell, debug etc..)
