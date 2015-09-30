from astropy.table import Table, Column
from scipy.stats.mstats import mode
import csv
import numpy as np
import matplotlib.pyplot as plt
import string
import pdb
import ast
import os




# using the NodeTransformer, you can also modify the nodes in the tree,
# however in this example NodeVisitor could do as we are raising exceptions
# only.
class Transformer(ast.NodeTransformer):
    ALLOWED_NAMES = set(['Decimal', 'None', 'False', 'True'])
    ALLOWED_NODE_TYPES = set([
        'Expression', # a top node for an expression
        'Tuple',      # makes a tuple
        'Call',       # a function call (hint, Decimal())
        'Name',       # an identifier...
        'Load',       # loads a value of a variable with given identifier
        'Str',        # a string literal

        'Num',        # allow numbers too
        'List',       # and list literals
        'Dict',       # and dicts...
    ])

    def visit_Name(self, node):
        if not node.id in self.ALLOWED_NAMES:
            raise RuntimeError("Name access to %s is not allowed" % node.id)

        # traverse to child nodes
        return self.generic_visit(node)

    def generic_visit(self, node):
        nodetype = type(node).__name__
        if nodetype not in self.ALLOWED_NODE_TYPES:
            raise RuntimeError("Invalid expression: %s not allowed" % nodetype)

        return ast.NodeTransformer.generic_visit(self, node)


def parse_panoptes(filename):
    user_name, user_ip = [], []
    workflow_id, created_at = [], []
    metadata, subject_imgname = [], []
    subject_id, subject_num = [], []
    value, workflowversion = [], []
    
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:

            # get rid of null values
            try:
                val = ast.literal_eval(row['annotations'])[0]['value'][0]
            except:
                val = ast.literal_eval(row['annotations'])[0]['value']
            if val != 'n':
                value.append(int(val))
                            
                user_name.append(row['user_name'])
                user_ip.append(row['user_ip'])
                workflow_id.append(row['workflow_id'])
                workflowversion.append(row['workflow_version'])
                created_at.append(row['created_at'])
                
                bigdict = ast.literal_eval(row['subject_data'])
                for k,v in bigdict.iteritems():
                    
                    subject_id.append(int(bigdict[k]['dr7objid']))
                    
                    try:
                        subject_imgname.append(bigdict[k]['imgname'])
                    except:
                        subject_imgname.append(bigdict[k]['imagename'])

                    subject_num.append(int(k))

    #pdb.set_trace()
    data = Table(data=(user_name, user_ip, workflow_id, created_at,  
                       subject_id, subject_num, subject_imgname, value, 
                       workflowversion),
                 names=('user_name', 'user_ip', 'workflow_id', 'created_at', 
                        'dr7objid', 'panoptes_id', 'imgname', 
                        'value', 'workflow_version'))
    name = os.path.splitext(filename)
    data.write(name[0]+'_output_parsed.txt',format='ascii')
    data.write(name[0]+'_output_parsed.fits', overwrite=True)
    return data

def fix_dr7objids(data, correctvals):
  
    for m in correctvals:
        m['imgname_2'] = m['imgname_2'].strip()
        
    for dat in data:
        match = np.where(dat['imgname']==matched['imgname_2'])
        #pdb.set_trace()
        
        dat['subject_id']=matched['dr7objid'][match]  

    return data

def explore_galaxy_classifications(data, unique_users, unique_gals, filename):

    number_classifications, number_null = [], []
    max_classification_single_classifier = []
    number_unique_classifiers = []
    median_all_vals, mean, std = [], [], []
    over_classified = []
    user_mode, median_user_mode, median_user_last = [], [], []
    final_val = []
    lucy_group, claudia_group = [], []

    for gal in unique_gals:
        mm = np.where(data['dr7objid']==gal)
        dats = data[mm]
        classifications = data['value'][mm]
        classifiers = data['user_name'][mm]

        # number of actual classifications
        number_classifications.append(len(classifications))
        
        median_all_vals.append(np.median(classifications))
        
        # how many unique classifiers for this galaxy? 
        unique_classifiers, counts = np.unique(classifiers, return_counts=True)
        
        # number of unique classifiers
        number_unique_classifiers.append(len(unique_classifiers))
   
        # number of classifications by the classifier 
        # with the most classifications
        max_classification_single_classifier.append(max(counts))
        
        # were those classifications all the same or different?
        class_vals_mode, class_vals_last = [], []
        lgroup, cgroup = [], []
        if len(unique_classifiers) > 1:
            for u in unique_classifiers:
                userclass_idx = np.where(classifiers == u)[0]
                userclass_vals = classifications[userclass_idx]
                
                # how do we choose which classification PER CLASSIFIER? 
                # There should only be one vote per person!
                (val_mode, count) = mode(userclass_vals, axis=None)
                class_vals_mode.append(val_mode[0])
                class_vals_last.append(userclass_vals[-1])

                if u.strip() in ['Mel23', 'KWillett', 'not logged in']:
                    lgroup.append(val_mode[0])
                if u.strip() in ['highzgal', 'aliwaffles', 
                                 'claudiascarlata2']:
                    cgroup.append(val_mode[0])

            #print class_vals_mode
            mode_vote, counts_vote = mode(class_vals_mode, axis=None)

            if (mode_vote == 1) and counts_vote > 1:
                final_val.append(1.0)
            elif (mode_vote == 1) and counts_vote == 1:
                final_val.append(np.median(class_vals_mode))
            elif np.median(class_vals_mode) == 1:
                final_val.append(5.0)
            else:
                final_val.append(np.median(class_vals_mode))

            #median_user_last.append(np.median(class_vals_last))
            lucy_group.append(np.median(lgroup))
            claudia_group.append(np.median(cgroup))

    final_val = np.array(final_val)
    # how many galaxies had only one classifier (namely, ME)
    more_than_one = np.where(np.array(number_unique_classifiers) != 1)
  
    class_total = Table(data=(unique_gals[more_than_one], final_val, 
                              lucy_group, claudia_group), 
                        names=('dr7objid', 'class', 'l_class', 'c_class'))
    class_total.write(filename+'_classifications_mto.fits', 
                      overwrite=True)

    pdb.set_trace()
    # compare galaxy classifications using ALL classifications by all users, 
    # the median value from each user for each galaxy, and the first answer 
    # given by each user for each galaxy
    
    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(121)
    ax1.hist(np.array(median_all_vals)[more_than_one])
    ax1.set_xlabel('Median of all classifications for each galaxy')
    '''
    ax2 = fig.add_subplot(122)
    ax2.hist(median_user_mode, color='green', alpha=0.5, label='Mode')
    ax2.hist(median_user_last, color='yellow', alpha=0.5, label='Last')
    ax2.set_xlabel('Median classification per user for each galaxy')
    ax2.legend(loc='best')
    plt.savefig(filename+'_classificaction_distribution_mto.png')
    #plt.show()
    plt.close()
    '''
    plt.figure()
    plt.hist(max_classification_single_classifier)
    plt.xlabel('Number of classifications made by the classifier with the most votes for a given galaxy')
    plt.savefig(filename+'_maxclassifier_mto.png')
    #plt.show()
    plt.close()

    plt.figure()
    plt.plot(number_classifications, number_unique_classifiers, 'r^')
    plt.xlabel('# classifications')
    plt.ylabel('# of unique classifiers')
    plt.savefig(filename+'_unique_classifiers_mto.png')
    #plt.show()
    plt.close()

    return classifications_final

def compare_classifications(data, filename, suffix):
    gzclass = data['zclass'].astype(int)
    panclass = data['class']
    l_class = data['l_class']
    c_class = data['c_class']
    
    classifications, l_classes, c_classes = [], [], []
    for p,l,c in zip(panclass, l_class, c_class):
        if (p <= 4.) & (p > 1.):
            classifications.append(4.0)
        if p == 1:
            classifications.append(1.0)
        if p < 1.0:
            classifications.append(2.0)
        if p == 5.0:
            classifications.append(3.0)

        if (l <= 4.) & (l > 1.):
            l_classes.append(4.0)
        if l == 1:
            l_classes.append(1.0)
        if l < 1.0:
            l_classes.append(2.0)
        if l == 5.0:
            l_classes.append(3.0)

        if (c <= 4.) & (c > 1.):
            c_classes.append(4.0)
        if c == 1:
            c_classes.append(1.0)
        if c < 1.0:
            c_classes.append(2.0)
        if c == 5.0:
            c_classes.append(3.0)
        #if (p <= 4.) & (p >= 2.):
        #    classifications.append(4.0)
        #if (p > 1.25) & (p < 2.):
        #    classifications.append(3.0)
        #if (p > .75) & (p <= 1.25):
        #    classifications.append(1.0)
        #if (p >= 0.) & (p <= 0.75):
        #    classifications.append(2.0)

    data['panclass'] = classifications
    #data['l_panclass']=l_classes
    #data['c_panclass']=c_classes
    data.write(filename+'_classifications_'+suffix+'_final_groupsplit.fits', 
               overwrite=True)
            
    fig = plt.figure()

    labels1 = ['mergers', 'ell', 'edge on', 'disks']
    labels2 = ['mergers', 'ell', '???', 'disks']

    ax1=fig.add_subplot(221)
    ax1.hist(gzclass, normed=1)
    ax1.set_ylim(0, 3)
    ax1.set_title('classifications from gz2')
    ax1.xaxis.set_ticks([1.0, 2.0, 3.0, 4.0])
    ax1.set_xticklabels(labels1)

    ax2=fig.add_subplot(222)
    ax2.hist(classifications, normed=1)
    ax2.set_ylim(0, 3)
    ax2.set_title('classifications from panoptes')
    ax2.xaxis.set_ticks([1.0, 2.0, 3.0, 4.0])
    ax2.set_xticklabels(labels2)

    ax3=fig.add_subplot(223)
    ax3.hist(l_classes, normed=1)
    ax3.set_ylim(0, 3)
    ax3.set_title("classifications from Lucy's group")
    ax3.xaxis.set_ticks([1.0, 2.0, 3.0, 4.0])
    ax3.set_xticklabels(labels2)

    ax4=fig.add_subplot(224)
    ax4.hist(c_classes, normed=1)
    ax4.set_ylim(0, 3)
    ax4.set_title("classifications from Claudia's group")
    ax4.xaxis.set_ticks([1.0, 2.0, 3.0, 4.0])
    ax4.set_xticklabels(labels2)

    plt.tight_layout()
    plt.savefig(filename+'_classcompare_gz_vs_panoptes_'+suffix+'_groupsplit.png')

    plt.show()
    plt.close()
    
    equal_count = 0
    #np.empty(())
    mergers, ellipticals, throwaway, disks = [], [], [], []
    for c1, c2 in zip(classifications, gzclass):
        if c1 == c2:
            equal_count+=1
            #equal.append(data['subject_id'])
        else:
            if c2 == 1.0:
                mergers.append(c1)
            if c2 == 2.0:
                ellipticals.append(c1)
            if c2 == 3.0:
                throwaway.append(c1)
            if c2 == 4.0:
                disks.append(c1)
        #pdb.set_trace()

    fig = plt.figure()
    ax=fig.add_subplot(111)
    ax.hist(classifications, normed=0, alpha=0.5, label='original gz2')
    ax.hist(ellipticals,color='orange', alpha=0.5, label='from ellipticals')
    ax.hist(throwaway,color='yellow',alpha=0.5, label="from 'edgeon'")
    ax.hist(disks,color='red',alpha=0.5, label='from disks')
    plt.savefig(filename+'_scatteredout_gz_panoptes_'+suffix+'.png')
    plt.legend(loc='upper left')
    plt.show()

#---------------------------------------------------------------------------#

def main():

    filename = 'panoptes_run2.csv'
    name = 'panoptes_run2'

    #try:
    #data = parse_panoptes(filename)
    #except:
    data = Table.read(name+'_output_parsed_matched.fits')
    
    workflows = np.unique(data['workflow_id'])
    '''
    for flow in workflows:
        thing = data[data['workflow_id']==flow]

        # get some overall stats for that run
        unique_users = np.unique(thing['user_name'])
        unique_gals = np.unique(thing['dr7objid'])
        
        print unique_users
        print len(unique_gals)
        #pdb.set_trace()
        classifications = explore_galaxy_classifications(thing, unique_users, 
                                                         unique_gals, name)
        #explore_user_classifications(unique_users, unique_gals)
        pdb.set_trace()
                                                         
    #'''

    data_class = Table.read(name+'_classifications_mto_final.fits')
    print len(data_class)

    #pdb.set_trace()
    compare_classifications(data_class, name, suffix='mto')


if __name__ == '__main__':
    main()
        
