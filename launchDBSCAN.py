"""
DBSCAN Launcher.
    :param 	dataFile: 	relative path to the file where the dataset is stored.
    :param	fragSize: 	percentage of the dynamic range for each partition.
    :param	minPoints:	minimum number of neighbours for a point to be considered a core point. 
    :param	*args: 		if the dataset is 2D or 3D, plotting can be enabled with a 2 or a 3 as an argument.
"""	
import sys
import DBSCAN_Adap2 as DBSCAN
import numpy as np

def main(dataFile, fragSize, epsilon, minPoints, numParts, *args):
    [defCluster, fragVec]=DBSCAN.DBSCAN(dataFile, int(fragSize), float(epsilon), int(minPoints), int(numParts))


    #Plotting if selected.
    newName=dataFile[dataFile.find('.')+1:dataFile.rfind('.')]
    newName='.'+newName+'n.txt'
    dataset=np.loadtxt(newName)
    if len(args)>0 and args[0]=='2D':
        #Try and plot regions division
        import matplotlib.pyplot as plt
        from matplotlib import colors	
        fig, ax = plt.subplots()
        plt.hlines(fragVec[1],0,1, 'k', 'dashdot', linewidths=0)
        plt.vlines(fragVec[0],0,1, 'k', 'dashdot', linewidths=0) 
        ax.scatter([p[0] for p in dataset], [p[1] for p in dataset], s=1)
        plt.savefig('dataset.png')
        plt.close()
        colours = [hex for (name, hex) in colors.cnames.iteritems()]
        fig, ax = plt.subplots()
        plt.hlines(fragVec[1],0,1, 'k', 'dashdot', linewidths=0.1)
        plt.vlines(fragVec[0],0,1, 'k', 'dashdot', linewidths=0.1) 
        for i,key in enumerate(defCluster):
                ax.scatter([p[0] for p in defCluster[i]], [p[1] for p in defCluster[i]], color=colours[i], s=1)
        plt.savefig('clusters.png')
        plt.close()
    elif len(args) > 0 and args[0]=='3D':
        #from mpl_toolkits.mplot3d import Axes3D
        fig=plt.figure()
        ax=fig.add_subplot(111, projection='3d')
        #fig, ax = plt.subplots(projection='3d')
        #ax.set_xticks(np.arange(0, 2, float(fragSize)))
        #ax.set_yticks(np.arange(0, 2, float(fragSize)))
        ax.scatter([p[0] for p in dataset], [p[1] for p in dataset], [p[2] for p in dataset], s=1)
        plt.grid()
        plt.savefig('plot/dataset.png')
        plt.close()
        colours = [hex for (name, hex) in colors.cnames.iteritems()]
        fig=plt.figure() 
        ax = fig.add_subplot(111, projection='3d')
        #ax.set_xticks(np.arange(0, 2, fragSize))
        #ax.set_yticks(np.arange(0, 2, fragSize))
        #ax.tick_params(colors='lightgrey')
        #For debuggin it might be usefull to let ticks on, otherwise they might get really annoying.
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        for c in defCluster:
                ax.scatter([p[0] for p in defCluster[c]], [p[1] for p in defCluster[c]], [p[2] for p in defCluster[c]], color=colours[c], s=1)
        #plt.grid()
        plt.show()
        #plt.savefig('plot/clusters.png')
        #plt.close()
    elif len(args) > 0 and args[0]=='text':
        #from sklearn.manifold import TSNE
        #import matplotlib.pyplot as plt
        #import gzip
        #import pickle

        # Restore the object
        with gzip.open('idWordVec.pklz', 'rb') as f:
                id_word_vec = pickle.load(f)

        tsne = TSNE(n_components=2)
        X_tsne = tsne.fit_transform([v for _, _, v in id_word_vec[0:2000]])
        plt.scatter(X_tsne[:, 0], X_tsne[:, 1])

        for i, word in enumerate([word for _, word, _ in id_word_vec[0:2000]]):
                plt.annotate(word, (X_tsne[i, 0], X_tsne[i, 1]))

        plt.show()
    f= open('outDBSCAN.txt', 'w')
    for num,lista in enumerate(defCluster):
        f.write('Cluster '+str(num)+':\n')
        for point in defCluster[num]:
                f.write(str(point)+'\n')	
    f.close()
    

if __name__=='__main__':
    main(sys.argv[1],float(sys.argv[2]),float(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6] if len(sys.argv) >= 7 else 0)
