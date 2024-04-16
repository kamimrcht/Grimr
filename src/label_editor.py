import csv
import pickle

def NOT(label: str, labels_dict: dict[str, set[int]], name: str = None):
    """
    Crée une nouvelle étiquette dans `labels_dict`, qui est le complémentaire de `label` dans [1,..,N]
    """
    if name is None:
        name = "NOT-" + label

    new_set = labels_dict["universe"].difference(labels_dict[label])
    if len(new_set) > 0:
        labels_dict[name] = new_set
    else:
        print(name, "is empty.")


def UNION(label_list: list[str], labels_dict: dict[str, set[int]], name: str = None):
    """
    Crée une nouvelle étiquette qui est l'union des étiquettes de `label_list`
    """
    if name is None:
        name = ""
        for label in label_list:
            name += label + "u"
        name = name[:-1]

    set_list = [labels_dict[label] for label in label_list]
    new_set = set.union(*set_list)
    if len(new_set) > 0:
        labels_dict[name] = new_set
    else:
        print(name, "is empty.")


def INTERSECTION(label_list: list[str], labels_dict: dict[str, set[int]], name: str = None):
    """
    Calcule l'intersection de `label_list` et crée une nouvelle étiquette si l'intersection est non vide
    """
    if name is None:
        name = ""
        for label in label_list:
            name += label + "n"
        name = name[:-1]

    set_list = [labels_dict[label] for label in label_list]
    new_set = set.intersection(*set_list)
    if len(new_set) > 0:
        labels_dict[name] = new_set
    else:
        print(name, "is empty.")


def DIFFERENCE(labelA: str, labelB: str, labels_dict: dict[str, set[int]], name: str = None):
    """
    Crée une nouvelle étiquette qui est la différence de `labelA` avec `labelB` : A\B
    """
    if name is None:
        name = labelA + "\\" + labelB

    new_set = labels_dict[labelA].difference(labels_dict[labelB])
    if len(new_set) > 0:
        labels_dict[name] = new_set
    else:
        print(name, "is empty.")


def SYMMETRIC_DIFFERENCE(labelA: str, labelB: str, labels_dict: dict[str, set[int]], name: str = None):
    """
    Crée une nouvelle étiquette qui est la différence symétrique de `labelA` et `labelB`
    """
    if name is None:
        name = labelA + "Δ" + labelB
    new_set = labels_dict[labelA].symmetric_difference(labels_dict[labelB])
    if len(new_set) > 0:
        labels_dict[name] = new_set
    else:
        print(name, "is empty.")


# LOAD METADATA

csv_file = '../test_files/metadata.csv'

labels_dict = dict()
map = dict()
universe=set()

with open(csv_file, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')

    index = 0

    for row in spamreader:

        filename = row[0]
        labels= row[1:]

        map[index]=filename
        universe.add(index)

        for label in labels:
            if label not in labels_dict.keys():
                labels_dict[label]=set()
            labels_dict[label].add(index)

        index+=1

    labels_dict['universe']=universe

# BUILD NEW LABELS

labels = [i for i in labels_dict.keys() if i != "universe"]

for label in labels:
    NOT(label, labels_dict)

UNION(['A', 'B', 'C'], labels_dict)
INTERSECTION(['A', 'G'], labels_dict)
DIFFERENCE("A", "H", labels_dict)
SYMMETRIC_DIFFERENCE("A", "H", labels_dict)

for label in labels_dict.keys():
    if label != "universe":
        print(label, labels_dict[label])

# SAVE LABELS

with open('../test_files/metadata.pickle', 'wb') as fichier:
    pickle.dump(labels_dict, fichier)

with open('../test_files/map.txt','a') as fichier:
    for key in map.keys():
        fichier.write(str(key)+'\t'+map[key]+'\n')