from astropy.table import Table
from scipy.stats.mstats import mode
import csv
import numpy as np
import matplotlib.pyplot as plt
import string
import pdb
import ast




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


def parse_panoptes():
    user_id, user_ip = [], []
    workflow_id, created_at = [], []
    metadata, subject_imgname = [], []
    subject_id, subject_num = [], []
    value, workflowversion = [], []
    
    counter=0
    with open('Panoptes_datadownload2.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            
            user_id.append(row['user_id'])
            user_ip.append(row['user_ip'])
            workflow_id.append(row['workflow_id'])
            workflowversion.append(row['workflow_version'])
            created_at.append(row['created_at'])
            
            bigdict = ast.literal_eval(row['subject_data'])
            for k,v in bigdict.iteritems():
                subject_id.append(bigdict[k]['# dr7objid'])
                subject_imgname.append(bigdict[k]['imgname'])
                subject_num.append(int(k))
                
            metadata.append(ast.literal_eval(row['metadata']))
            
            ##try:
            #value.append(ast.literal_eval(row['annotations'])[0]['value'][0])
            #except:
            value.append(ast.literal_eval(row['annotations'])[0]['value'])

            counter+=1

        
    data = Table(data=(user_id, user_ip, workflow_id, created_at,  
                       subject_id, subject_num, subject_imgname, value, 
                       workflowversion),
                 names=('user_id', 'user_ip', 'workflow_id', 'created_at', 
                        'subject_id', 'subject_num', 'imgname', 
                        'value', 'workflow_version'))
    data.write('Panoptes_output_parsed.txt',format='ascii')


def fix_dr7objids(data, correctvals):
  
    for m in correctvals:
        m['imgname_2'] = m['imgname_2'].strip()
        
    for dat in data:
        match = np.where(dat['imgname']==matched['imgname_2'])
        #pdb.set_trace()
        
        dat['subject_id']=matched['# dr7objid'][match]  

    return data

def explore_galaxy_classifications(unique_users, unique_gals):
    number_classifications, number_null = [], []
    max_classification_single_classifier = []
    number_unique_classifiers = []
    median_all_vals, mean, std = [], [], []
    over_classified = []
    median_user_mode, median_user_first = [], []

    for gal in unique_gals:
        mm = np.where(thing['subject_id']==gal)
        dats = thing[mm]
        classification_vals = thing['value'][mm]
        classifiers = thing['user_id'][mm]

        # Are there "null" classificaction values?
        #if np.any(classification_vals == 'n'):
        notnulls = np.where(classification_vals != 'n')
        
        # number of "null" classifications
        number_null.append(len(classification_vals) - len(notnulls[0]))

        # number of actual classifications
        number_classifications.append(len(notnulls[0]))
        
        # values of the actual classifications
        class_vals = np.array([int(c) for c in \
                               classification_vals[notnulls]])
        median_all_vals.append(np.median(class_vals))
        
        # how many unique classifiers for this galaxy? 
        uni, counts = np.unique(classifiers[notnulls], return_counts=True)
        
        # number of unique classifiers
        number_unique_classifiers.append(len(uni))
        
        # number of classifications by the classifier 
        # with the most classifications
        max_classification_single_classifier.append(max(counts))
        
        # were those classifications all the same or different?
        class_vals_mode, class_vals_first = [], []
        for u in uni:
            userclass_idx = np.where(classifiers[notnulls] == u)[0]
            userclass_vals = class_vals[userclass_idx]

            # how do we choose which classification PER CLASSIFIER? 
            # There should only be one vote per person!
            (val_mode, count) = mode(userclass_vals, axis=None)
            class_vals_mode.append(val_mode[0])
            class_vals_first.append(userclass_vals[0])

        median_user_mode.append(np.median(class_vals_mode))
        median_user_first.append(np.median(class_vals_first))
        
        #if len(classification_vals) > 20:
            #for u, v in zip(classifiers, classification_vals):
            #    print u, v
            #pdb.set_trace()

    pdb.set_trace()
    classifications_final = Table(data=(unique_gals, median_user_mode), 
                                  names=('gals', 'class'))
    classifications_final.write('panoptes_classifications_final.fits', 
                                overwrite=True)
    pdb.set_trace()

    # compare galaxy classifications using ALL classifications by all users, 
    # the median value from each user for each galaxy, and the first answer 
    # given by each user for each galaxy
    
    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(121)
    ax1.hist(median_all_vals)
    ax1.set_xlabel('Median of all classifications for each galaxy')

    ax2 = fig.add_subplot(122)
    ax2.hist(median_user_mode, color='green', alpha=0.5, label='Mode')
    ax2.hist(median_user_first, color='yellow', alpha=0.5, label='First')
    ax2.set_xlabel('Median of the XXX classifications per user for each galaxy')
    ax2.legend(loc='best')
    plt.savefig('panoptes_galaxy_classificactions.png')
    #plt.show()
    plt.close()

    plt.figure()
    plt.hist(max_classification_single_classifier)
    plt.xlabel('Number of classifications made by the classifier with the most votes for a given galaxy')
    plt.savefig('panoptes_max_classifier.png')
    #plt.show()
    plt.close()

    plt.figure()
    plt.plot(number_classifications, number_unique_classifiers, 'r^')
    plt.xlabel('# classifications')
    plt.ylabel('# of unique classifiers')
    plt.savefig('panoptes_unique_classifiers.png')
    #plt.show()
    plt.close()

    plt.figure()
    plt.plot(number_classifications, number_null, 'ko')
    plt.xlabel('# classifications')
    plt.ylabel("# 'null' classifications")
    plt.savefig('panoptes_null_total_classification.png')
    #plt.show()
    plt.close()


#---------------------------------------------------------------------------#

data = Table.read('panoptes_gz2_classification_match.fits')

gzclass = data['zclass'].astype(int)
#gzclass = data[]
panclass = data['class']

classifications=[]
for p in panclass:
    if (p < 4.) & (p >= 2.):
        classifications.append(4.0)
    if (p > 1.) & (p < 2.):
        classifications.append(3.0)
    if (p > .5) & (p <= 1.):
        classifications.append(1.0)
    if (p >= 0.) & (p <= 0.5):
        classifications.append(2.0)

fig = plt.figure()
ax1=fig.add_subplot(121)
ax1.hist(gzclass)
ax1.set_ylim(0, 300)
ax1.set_title('classifications from gz2')
ax2=fig.add_subplot(122)
ax2.hist(classifications)
ax2.set_ylim(0, 300)
ax2.set_title('classifications from US')
plt.savefig('panoptes_gz_class_compare.png')
#plt.show()
plt.close()

equal_count = 0
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
ax.hist(gzclass,alpha=0.5, label='original gz2')
ax.hist(ellipticals,color='orange', alpha=0.5, label='from ellipticals')
ax.hist(throwaway,color='yellow',alpha=0.5, label="from 'edgeon'")
ax.hist(disks,color='red',alpha=0.5, label='from disks')
plt.savefig('panoptes_gzclass_scatteredout.png')
plt.legend(loc='upper left')
plt.show()

exit()

#parse_panoptes()

data = Table.read('Panoptes_run1_old/Panoptes_output_parsed.fits')
matched = Table.read('Panoptes_run1_old/Panoptes_output_dr7matched.fits')

# fix the objectIDs
data_fix = fix_dr7objids(data, matched)

data_fix.write('Panoptes_output_parsed.fits', overwrite=True)

workflows = np.unique(data['workflow_id'])

run1 = data[data['workflow_id']==workflows[0]]
run2 = data[data['workflow_id']==workflows[1]]

for thing in [run2]:
    # get some overall stats for that run
    unique_users = np.unique(thing['user_id'])
    unique_gals = np.unique(thing['subject_id'])

    explore_galaxy_classifications(unique_users, unique_gals)

    #explore_user_classifications(unique_users, unique_gals)

    for u in unique_users:
        mm = np.where(thing['user_id']==u)
        userclass
        
        userclass_idx = np.where(classifiers[notnulls] == u)[0]
        userclass_vals = class_vals[userclass_idx]
        
        # how do we choose which classification PER CLASSIFIER? 
        # There should only be one vote per person!
        for gal in unique_gals:
            (val_mode, count) = mode(userclass_vals, axis=None)
            class_vals_mode.append(val_mode[0])
            class_vals_first.append(userclass_vals[0])
            
            median_user_mode.append(np.median(class_vals_mode))
            median_user_first.append(np.median(class_vals_first))
            
            
        classification_vals = thing['value'][mm]
        classifiers = thing['user_id'][mm]
        
        
 
pdb.set_trace()

        
