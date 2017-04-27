#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Thibault TUBIANA"
__version__  = "3.9.1"
__copyright__ = "copyleft"
__license__ = "GNU GPLv3"
__date__ = "2016/11"

#==============================================================================
#                     MODULES                            
#==============================================================================
import argparse
import mdtraj as md
import numpy as np
import os
import sys
import progressbar as pg
import datetime
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.cluster.hierarchy as sch
import operator
from prettytable import PrettyTable

try:
    import argcomplete
except:
    print("argcomplete not detected, if you want to use autocompletetion")
    print("install argcomplete module.")
    print("See https://github.com/kislyuk/argcomplete for more info")
    pass
#==============================================================================
#                     GLOBAL VARIABLES                            
#==============================================================================
WIDGETS = [pg.Bar('>'), ' ', pg.ETA(), ' ', pg.ReverseBar('<')]
COORDS=[]
COLOR_LIST = ["red","blue","lime","gold",
              "darkorchid", "deepskyblue",
              "orange","brown", "gray","black",
              "darkgreen","navy"]

class Cluster_class():
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

#==============================================================================
#                     TOOL FONCTIONS
#==============================================================================
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
    LOGFILE.write("command line     : python ")
    i = 0
    while i < len(sys.argv):
        LOGFILE.write("{} ".format(sys.argv[i]))
        if sys.argv[i] in ["-st","-sr","-sa"]:
            i+=1
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
    
    LOGFILE.write("========================================================\n")
    LOGFILE.write("============  3D STRUCTURES CLUSTERING {}  ============\n"\
        .format(__version__))
    LOGFILE.write("========================================================\n")
    LOGFILE.write("\n")
    
    
    LOGFILE.write("************ General information ************\n")
    LOGFILE.write("software version : {}\n".format(__version__))
    LOGFILE.write("Created on       : {}\n".format(datetime.datetime.now()))
    write_command_line()
    LOGFILE.write("ARGUMENTS : \n")
    LOGFILE.write("  Selection string :\n")
    LOGFILE.write("      Atoms selected in trajectory = {} \n".format(
                                                        selection_string))
    LOGFILE.write("      Atoms selected for alignement = {} \n".format(
                                                        select_align))
    LOGFILE.write("      Atoms selected for RMSD = {} \n".format(select_rmsd))
    LOGFILE.write("  trajectory file  : {} \n".format(traj))
    LOGFILE.write("   Number of frames  : {} \n".format(mdtrajectory.n_frames))
    LOGFILE.write("   Number of atoms  : {} \n".format(mdtrajectory.n_atoms))
    LOGFILE.write("  topology file    : {} \n".format(topo))
    LOGFILE.write("  method used of clusterring : {}".format(args["method"]))
    LOGFILE.write("\n\n")
    if args["ngroup"]:
        LOGFILE.write("  Number of cluster asked: {}\n".format(args["ngroup"]))
    if args["cutoff"]:
        LOGFILE.write("  cutoff for dendrogram clustering: {}\n".format("cutoff"))
        
        
def extract_selected_atoms(selection, traj):
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
    Returns:
        subtraj (mdtraj.trajecotry): subtrajectory of selected atoms
    """
    try:
        subtraj=traj.atom_slice(traj.top.select(selection))
        subtraj.center_coordinates()
        return subtraj
    except:
        print("ERROR : there is an error with your selection string")
        print("        SELECTION STRING : ")
        print(("        {}".format(selection)))
        print("        > Please check 'http://mdtraj.o\
                        rg/latest/atom_selection.html'")
        exit(1)
    
def send_error_message(calc_type, selection_string):
    """
    DESCRIPTION
    Print information regarding your selection string if this one is not 
    recognized by mdtraj
    Args:
        calc_type
    """
    print(("ERROR : {} selection string not valid".format(calc_type)))
    print(("        >{}".format(selection_string)))
    exit(1)   

def return_selection_atom(use_for,traj, selection_string):
    """
    DESCRIPTION
    return indices of selected atoms.
    ----
    Args:
        use_for (string): witch selection string was wrong ? (sr/sa/sr)
        traj (mdtraj.trajectory): trajectory
        selection_string (string): selection string wich produce an error
    """
    try:
        selection=traj.top.select(selection_string)
    except:
        send_error_message(use_for,selection_string)
    
    if len(selection)==0:
        send_error_message(use_for,selection_string)
    else:
        return selection
        
def save_dist_mat(distmat, rmsd_string):
    """
    DESCRIPTION
    Save the numpy matrix to reused afterward
    ----
    Args:
        distmat (numpy matrix): distance matrix 
        alignement (str) : alignement string, used for the name
    return:
        None
    """
    if rmsd_string:
        name=rmsd_string.replace(" ","_")
    else:
        name="matrix_all"
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
    returns: 
        clusters_labels_order(list): list of reorder clusters label
    """
    dict_order = {}
    for cluster in clusters:
        dict_order[cluster.id] = cluster.frames[0]
    #Evaluate order
    sorted_clusters = sorted(dict_order.items(), key=operator.itemgetter(1))
    for i in range(len(sorted_clusters)):
        dict_order[sorted_clusters[i][0]] = i+1 #  i+1 == reorder cluster number
    #reordering
    for cluster in clusters:
        cluster.id = dict_order[cluster.id]
                
#==============================================================================
#                     FONCTIONS
#==============================================================================


def parseArg():
    """
    This fonction will the list of pdb files and the distance
    @return: dictionnary of arguments
    Ex : 
    python Cluster_Analysis.py -f *.pdb -s A:1-30:CA
    """
    arguments=argparse.ArgumentParser(description="\
          This program was developped in order to clusterize molecular dynamic\
          trajectories. Amber, gromacs, chamm, namd, PDB")
    try:
        argcomplete.autocomplete(arguments)
    except:
        print("Warning : argcomplete module not detected")
        pass
    arguments.add_argument('-f', "--traj", help="trajectory file", required=True)
    arguments.add_argument('-t','--top', help="topfile", default=None)
    arguments.add_argument('-l','--logfile', help="logfile (logfile.txt). The \
        name of your output file will be the basename (name before the extention\
        of this logfile", default="logfile.txt")
    arguments.add_argument('-st','--select_traj', help="selection syntaxe for\
        Don't forget to add QUOTES besite this selection string.\
        trajectory extraction (all).", default="all")
    arguments.add_argument('-sa','--select_alignement', help="selection syntaxe\
        for alignement (backbone). Don't forget to add QUOTES besite this \
        selection string.", default="backbone")
    arguments.add_argument('-sr','--select_rmsd', help="selection syntaxe for \
    RMSD. Don't forget to add QUOTES besite this selection string.", default=None)
    
    #Clustering arguments
    arguments.add_argument('-m','--method', help="method for clustering : single\
       ; complete; average; weighted; centroid; median. (ward)", default="ward")
    arguments.add_argument('-cc',"--cutoff", help="cutoff for clusterization from\
                            hierarchical clusturing with Scipy", default=None)
    arguments.add_argument('-ng',"--ngroup", help="number of group asked. Use the maxclust method to clusterize in this case", default=None)                         

  
    #Interactive mode for distance matrix:
    arguments.add_argument('-i','--interactive', help="Interactive mode for distance matrix (Y/n)", default="Y")
    args = vars(arguments.parse_args())
    return(args)
    


def ask_choice(args,name):
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
        print("Interactive mode desactived. I will use the existing distance matrix")
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
        return(name)
    elif choice.upper() == "N":  # don't want to use it.. Recalculate!
        print("Calculation mode activated")
        return None
    elif choice.upper() == "O":  # I want to use another npy distance matrix
        npy_files=glob.glob("*.npy")
        for i,file in enumerate(npy_files):
            print("  {0} - {1}".format(i+1, file))
        print(" -->Please chooce and press Enter")
        #Check if the user give a good answer
        choice_file=input()
        try:
            name=npy_files[int(choice_file)-1]
            return name
        except:
            print("I didn't understand. Please try again")
            print("........")
            return ask_choice(args,name)               
    else:
        print("I didn't understand. Please try again")
        print("........")
        return ask_choice(args,name)

def search_dist_mat(rmsd_string, args):
    """
    Search if the distance matrix already exist
    ----
    Args:
        rmsd_string (str) : name of the numpy matrix
    """
    if rmsd_string:
        name=rmsd_string.replace(" ","_")
    else:
        name="matrix_all.npy"
    #Searching all npy file in the folder
    npy_files=glob.glob("*.npy")
    
    name = name+".npy"
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
        Clusters (): Clusters list
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
                #We don't want to calcul the same frame.
                if not frame_j == frame_i:
                    # we add to the corresponding value in the list of all rmsd
                    # the RMSD betwween frame_i and frame_j
                    mean_rmsd_per_frame[frame_i] += DM[frame_i-1,frame_j-1]
            # mean calculation
            mean_rmsd_per_frame[frame_i] /= len(frames)
            
            # Representative frame = frame with lower RMSD between all other
            # frame of the cluster
            repre = min(mean_rmsd_per_frame, key=mean_rmsd_per_frame.get)
            cluster.representative = repre+1 # Don't forget +1 to get the 
                                             # real frame number
            
            # spread = mean rmsd in all the cluster (*10 to have angstöm)
            cluster.spread = sum(mean_rmsd_per_frame.values()) / len(frames) 
            cluster.spread *= 10
        
        
def create_DM(traj, alignement_string, rmsd_string,args):
    """
    DESCRIPTION 
    Calcul the distance matrix
    ---
    Args:
        traj (mdtraj.trajectory): trajectory
        alignement_string (string): string for trajectory alignement
        rmsd_string (string): atom selection for rmsd calculation (and 
                              matrix distance calculation)
        args (dict): all arguments in dictionary
    return:
        distances (numpy matrix): distance matrix
    """
    #Get Atoms indices from selection string
    if rmsd_string:
        print("NOTE : Extraction of subtrajectory for time optimisation")
        traj = extract_selected_atoms(rmsd_string, traj)

    alignement_selection = return_selection_atom(use_for = "ALIGNEMENT",\
                                            traj   = traj,\
                                            selection_string=alignement_string)  
    
    # Trajectory superposition  (aligment)
    traj_aligned = traj.superpose(traj[0],
                                  atom_indices=alignement_selection,\
                                  parallel=True)
    
    # matrix initialization                              
    distances = np.empty((traj.n_frames, traj.n_frames))
    
    # Searching if a distance file already exist    
    distance_file=search_dist_mat(rmsd_string,args)
    
    # If a distance matrix file was found and choosed, we load it.
    if distance_file:
        return np.load(distance_file)
        printScreenLogfile(" >Distance Matrix File Loaded!")
    else:  # otherwise 
        pbar = pg.ProgressBar(widgets=WIDGETS, maxval=traj.n_frames).start()
        counter=0
        # Pairwise RMSD calculation (matrix n²) 
        for i in range(traj.n_frames):
            distances[i]=md.rmsd(traj_aligned,traj_aligned, frame=i)
            pbar.update(counter)
            counter+=1
        pbar.finish()
        
        #Finaly, we save the matrix if we want to load it again afterward
        print("Calculation ended - saving distance matrix")
        save_dist_mat(distances, rmsd_string)
        return distances       

def onclick(event):
    """
    DESCRIPTION
    This function is used to get coordinate of the mouse on the matplotlib
    windows.
    """
    ix, iy = event.xdata, event.ydata
    
    global COORDS
    COORDS.append((ix,iy))
    
    #juste one clic
    if len(COORDS)==1:
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
    #cluster_with_frame_number = defaultdict(lambda : [])
    clusters_list = []
    for cluster_num in set(labels):
        clusters_list.append(Cluster_class(cluster_num))  # create new instance of cluster
        
    for i, cluster_num in enumerate(labels):
        clusters_list[cluster_num-1].frames.append(i)
        # for DEBUG
        if cluster_num != clusters_list[cluster_num-1].id :
            print ("{0} - {0}".format(cluster_num, clusters_list[cluster_num-1]))
            sys.exit(1)
        
    for cluster in clusters_list:
        cluster.size = len(cluster.frames)
    #print(mapping)
    return clusters_list
    
    
    
def create_cluster_table(traj,args):
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
    select_align=args["select_alignement"]
    select_rmsd=args["select_rmsd"]
    cutoff=args["cutoff"]
    ncluster = args["ngroup"]
    output_graph_name = args["logfile"][:-4]
    #Creation of the distance matrix
    print("         creating distance matrix")
    distances=create_DM(traj, select_align, select_rmsd,args)
    if select_rmsd == None:
        select_rmsd = "None"
        
    linkage_file=search_dist_mat(select_rmsd+" linkage "+args["method"],args)
    
 
    if linkage_file:
        linkage = np.load(linkage_file)
    else:
        print("         Matrix shape: {}".format(distances.shape))
        print("         Scipy linkage in progress. Please wait. It can be long")
        print("         (approximatly 2mn30 for a 5000,5000 sized matrix)")
        try:
            linkage=sch.linkage(distances, method=args["method"])
        except:
            printScreenLogfile("ERROR : method name given for clustering didn't recognized")
            printScreenLogfile("      : methods are : single; complete; average; weighted; centroid; ward.")
            printScreenLogfile("      : check https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/sc"
            "ipy.cluster.hierarchy.linkage.html for more info")
            sys.exit(1)
        print("         >Done!")
        print("         ...Saving linkage matrix...")
        save_dist_mat(linkage,select_rmsd+" linkage "+args["method"])
        print("         >Done!")
        

    #if a cuttof for distance cuting is given
    if cutoff:
        clustering_result = sch.fcluster(linkage, cutoff, "distance")
    #otherwise we choose it on the screen by cliking on the matplotlib windows
    #If a number of wanted cluster is given
    elif ncluster:
        clustering_result = sch.fcluster(linkage,t=ncluster, criterion="maxclust")
        n_group=len(np.unique(clustering_result))
        cutoff = linkage[-(n_group-1),2]
    else:
        fig = plt.figure()
        fig.canvas.mpl_connect('button_press_event',onclick)
        plt.title("Please click where you wan to build clusters")
        sch.dendrogram(linkage)
        plt.show()
        cutoff=COORDS[0][1]
        clustering_result = sch.fcluster(linkage, cutoff, "distance")
    
    #write graphic
    fig = plt.figure()
     
    den = sch.dendrogram(linkage, color_threshold=float(cutoff))
    plt.axhline(y=int(cutoff), color = "grey")
    
    #Graph parameters
    plt.title("Clustering Dendrogram")
    ax = plt.axes() 
    ax.set_xticklabels([])
    ax.set_ylabel("Distance (AU)")
    ax.set_xlabel("Frames")
    
    plt.savefig("{}.png".format(output_graph_name), format="png", dpi=300, transparent=True)
        
    printScreenLogfile("  cutoff for clustering : {0}".format(cutoff))
    return distances,clustering_result


def write_representative_frame(traj, cluster):
    """
    DESCRIPTION
    Write representative frame of a cluster
    ----
    Args:
        traj (mdtraj.trajectory): trajectory
        cluster (Cluster_class): a Cluster object
    """
    if not os.path.exists("Cluster_PDB"):
        os.makedirs("Cluster_PDB")
        
    cluster_num = cluster.id
    frame = cluster.representative
    size = cluster.size
    traj[frame].save_pdb("Cluster_PDB/C%i-f%i-s%i.pdb" %(cluster_num,\
                                                                frame, size))

def create_linear_cluster_mapping_graph(clusters_list, output, size):
    """
    DESCRIPTION
    Create a linear cluster mapping graph where every frame is printed as a 
    colored barplot
    Args:
        clusters_labels (list): list of cluster number per frame
        output (string) output name for graph
    """
    global COLOR_LIST
    # order clusters_labels by order of appearance in the trajectory
    clusters_number_ordered = [0] * size
    for cluster in clusters_list:
        for frame in cluster.frames:
            clusters_number_ordered[frame] = cluster.id
    # DEFINE COLOR MAP
    # if two much cluster number to define colours by hand
    n_clusters = len(set(clusters_number_ordered))
    if n_clusters > 12 : 
        cmap = "hsv"
    else:
        # imshow take the last color for the last group (if 3 cluster, color of
        # clusters 3 will be brown")
        COLOR_LIST = COLOR_LIST[:n_clusters]
        cmap = mpl.colors.ListedColormap(COLOR_LIST)

    data = np.asmatrix(clusters_number_ordered)
    fig = plt.figure(figsize=(10,1))
    # move the graphic into the corner
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    # remove axes
    ax.set_axis_off()
    # set axes
    fig.add_axes(ax)
    # create graphic
    print(data.shape)
    ax.imshow(data,aspect='auto', interpolation='none', cmap=cmap)
    plt.savefig("{}-linear.png".format(output[:-4]), dpi=300)
    plt.close()

def get_RMSD_cross_cluster(clusters_list, distances):
    """
    DESCRIPTION
    This function will get the RMSD between all representativ frames of all
    clusters. Print it to the console and write it on the logfile
    Args:
        Clusters_list (list) : list of all clusters
        distances (np matrix) : rmsd matrix
    returns:
        None
    """
    # 1 - table preparation
    table = PrettyTable()
    n_clusters = len(clusters_list)
    field_names = ["Clusters"] + ["C"+str(x) for x in range(1,n_clusters+1)]
    table.field_names = field_names
    table.float_format = ".2"
    # + variable which countains all value except 0.0 values (for 
    # average calculation)
    non_diag_value = []
    
    # 2 - RMSD "calculation"
    RMSD_matrix = np.zeros((n_clusters,n_clusters))
    for i in range(n_clusters):
        repr1 = clusters_list[i].representative
        for j in range(n_clusters):
            repr2 = clusters_list[j].representative
            if i == j :
                RMSD_matrix[i][j] = 0.00
            else:
                #reprx-1 to correspond with the numpy index
                rmsd = distances[repr1-1][repr2-1]*10
                RMSD_matrix[i][j] = rmsd
                non_diag_value.append(rmsd)
                           
    
    # 3 - PrettyTable creation
    for i in range(n_clusters):
        table.add_row(["C"+str(i+1)] + RMSD_matrix[i].tolist())
        
    # 4 - print table
    printScreenLogfile("----------------------------")
    printScreenLogfile("RMSD MATRIX BETWEEN CLUSTERS")
    printScreenLogfile(table)
    printScreenLogfile("\nAVERAGE RSMD BETWEEN CLUSTERS : {:.2f}".format(
                                                np.mean(non_diag_value)))
    np.savetxt("RMSD_between_clusters.csv",RMSD_matrix,delimiter=";")
    

def Cluster_analysis_call(args):
    """
    DESCRIPTION
    Main function of the program : call other function as a pipeline
    Args:
        args (dict): all arguments in a dictionary
    """
    trajfile=args["traj"]    
    topfile=args["top"]
    select_traj=args["select_traj"]

    
    print("======= TRAJECTORY READING =======")
    if topfile == None and trajfile[-4:] == ".pdb":
        traj=md.load_pdb(trajfile)
    else:
        traj=md.load(trajfile,\
                     top=topfile)
                     
    init_log(args, traj)
    
    if not select_traj == "all":
        print("======= EXTRACTION OF SELECTED ATOMS =======")
        traj=extract_selected_atoms(select_traj, traj)      
    
    
    print("====== Clustering ========")
    distances,clusters_labels=create_cluster_table(traj,args)
    
    printScreenLogfile( "\n**** Cluster Results")
    clusters_list = return_mapping_cluster(clusters_labels)

    print("====== Reordering clusters ======")    
    reorder_cluster(clusters_list)
    # reordering the list by the cluster number
    clusters_list.sort(key = operator.attrgetter("id"))
    
    create_linear_cluster_mapping_graph(clusters_list, args["logfile"], traj.n_frames)
    calculate_representative_frame_spread(clusters_list, distances)
    
    for cluster in clusters_list:
          printScreenLogfile( "cluster {}".format(cluster.id))
          printScreenLogfile( "    size = {}".format(cluster.size))
          printScreenLogfile( "    representative frame={}".format(
            cluster.representative))
          printScreenLogfile( "    Frames : {} ".format(str([x+1 for x in cluster.frames])))
          printScreenLogfile( "    spread  : {} ".format(cluster.spread))
          write_representative_frame(traj, cluster)
  
    get_RMSD_cross_cluster(clusters_list, distances)



    
    
        
###############################################################################
#####                               MAIN                                 ######
###############################################################################
if __name__ == "__main__":
    print("********************************************************")
    print("**********  3D STRUCTURES CLUSTERING {} **************".format(\
              __version__))
    print("********************************************************")
    print("")
    #We get all arguments
    args=parseArg()
    global LOGFILE  
    LOGFILE=open("{}".format(args["logfile"]),"w")
    Cluster_analysis_call(args)
    LOGFILE.close()