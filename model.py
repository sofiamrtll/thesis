import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
import sklearn.metrics as metrics
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow import keras
import keras_tuner
from Bio import SeqIO
from pybedtools import BedTool
import matplotlib.pyplot as plt
import pylab as pl
import seaborn as sns


# One hot encoding
def fasta_to_onehotencode(seq):
    base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    sequence = seq
    sequence_int = [base2int.get(base, 9999) for base in sequence]  # Encode sequence bases as integers, i.e. A as 0, C as 1, etc.
    sequence_onehot = tf.one_hot(sequence_int, depth=4)
    return sequence_onehot


def fastatoarray(fasta_sequences):
    seq_array = np.zeros((1, 400, 4))
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        new_sequence = fasta_to_onehotencode(sequence)  # onehotencode(sequence)
        new_sequence = np.expand_dims(new_sequence, axis=0)
        seq_array = np.vstack((seq_array, new_sequence))
    seq_array = np.delete(seq_array, 0, 0)  # to remove the first array of zeros
    return seq_array


def get_df_m6A(path, file):
    m6A = BedTool(path+'/'+file)
    df_m6A = m6A.to_dataframe(disable_auto_names=True, names=[i for i in range(19)])
    df_m6A['abs_pos'] = df_m6A[1] - df_m6A[14]  # getting the absolute position of the m6A site in the sequence
    return df_m6A


def m6A_onehot(path, file, df_m6A):
    onehot_m6A_seq = np.zeros((1, 400, 5))
    m6A_file = np.zeros((1, 400, 1))
    seq_array = np.zeros((1, 400, 4))
    fasta_sequences = SeqIO.parse(open(path+'/'+file), 'fasta')
    #print('df_m6A')
    #print(df_m6A)
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        chrom = str(name).split(':')[0]
        start = (str(name).split(':')[1]).split('-')[0]
        end = ((str(name).split(':')[1]).split('-')[1]).split('(')[0]

        # filtering for the sequences which contain an m6A site
        filtered = df_m6A.loc[((df_m6A[13] == chrom) & (df_m6A[14] == int(start)) & (df_m6A[15] == int(end)))] #this could not be right 
        #print(filtered)
        # one hot encoding the sequence
        new_sequence = fasta_to_onehotencode(sequence)
        new_sequence = np.expand_dims(new_sequence, axis=0)
        seq_array = np.vstack((seq_array, new_sequence))
        m6A_seq = np.zeros((1, 400, 1))  # standard m6A column
        if filtered.empty is False:
            for i in filtered.index:
                index = filtered.loc[i, 'abs_pos'] - 1
                m6A_seq[0][index] = 1  # substituting with 1 in the m6A site
        m6A_file = np.vstack((m6A_file, m6A_seq))
    m6A_file = np.delete(m6A_file, 0, 0)
    seq_array = np.delete(seq_array, 0, 0)
    onehot_m6A_seq = np.c_[seq_array, m6A_file]
    return onehot_m6A_seq

def get_df_m6A_aug(path, file): 
    m6A = BedTool(path+'/'+file)
    df_m6A = m6A.to_dataframe(disable_auto_names= True, names = [i for i in range(36)])
    #print(df_m6A)
    df_m6A['abs_pos'] = df_m6A[1] - df_m6A[14]  # getting the absolute position of the m6A site in the sequence
    return df_m6A

def getonehotencoded_files(path, augmentation: bool): 
    '''
    augmentation = [True, False]
    '''
    directory = os.listdir(path)

    positive_m6A = np.zeros((1, 400, 5))
    neg1_m6A = np.zeros((1, 400, 5))
    neg2_m6A = np.zeros((1, 400, 5))

    df_pos, df_neg1, df_neg2 = get_miclip_filtered_df(path)

    if augmentation is True:
        df_pos_aug, df_neg1_aug, df_neg2_aug = get_miclip_filtered_df_aug(path)

    for file in directory:
        if 'fasta' in file and ('out' not in file):
            if ('positive' in file):
                positive_m6A = np.vstack((positive_m6A, m6A_onehot(path, file,
                                                                   df_pos)))

            if ('negative-1' in file):
                neg1_m6A = np.vstack((neg1_m6A, m6A_onehot(path, file,
                                                           df_neg1)))

            if ('negative-2' in file):
                neg2_m6A = np.vstack((neg2_m6A, m6A_onehot(path, file,
                                                           df_neg2)))
        if augmentation is True:
            if 'fasta' in file and ('out' in file):
                if ('positive' in file):
                    positive_m6A = np.vstack((positive_m6A, m6A_onehot(path, file,
                                                                       df_pos_aug)))
                if ('negative-1' in file):
                    neg1_m6A = np.vstack((neg1_m6A, m6A_onehot(path, file,
                                                               df_neg1_aug)))
                if ('negative-2' in file):
                    neg2_m6A = np.vstack((neg2_m6A, m6A_onehot(path, file,
                                                               df_neg2_aug)))
    positive_m6A = np.delete(positive_m6A, 0, 0)
    neg1_m6A = np.delete(neg1_m6A, 0, 0)
    neg2_m6A = np.delete(neg2_m6A, 0, 0)

    return positive_m6A, neg1_m6A, neg2_m6A  # for one fold


def get_miclip_filtered_df_aug(path):
    for file in os.listdir(path):
        if 'slop.miclip.filt.bed.out' in file:
            if ('positive' in file):
                df_pos_aug = get_df_m6A_aug(path, file)

            if ("negative-1" in file):
                df_neg1_aug = get_df_m6A_aug(path, file)

            if ("negative-2" in file):
                df_neg2_aug = get_df_m6A_aug(path, file)
    return df_pos_aug, df_neg1_aug, df_neg2_aug


def get_miclip_filtered_df(path):
    for file in os.listdir(path):  # all files
        if '.miclip.filt.bed.out' in file and ('slop' not in file):
            if ('positive' in file):
                df_pos = get_df_m6A(path, file)
            if ('negative-1' in file):
                df_neg1 = get_df_m6A(path, file)
            if ('negative-2' in file):
                df_neg2 = get_df_m6A(path, file)
    return df_pos, df_neg1, df_neg2


def getonehotencoded_files_pos(path):
    directory = os.listdir(path)

    df_pos = pd.DataFrame()
    df_pos_aug = pd.DataFrame()
    positive_m6A = np.zeros((1, 400, 5))

    for file in directory:
        if '.miclip.filt.bed.out' in file:
            if ('positive' in file) and ('slop' not in file):
                df_pos = get_df_m6A(path, file)
            if ('positive' in file) and ('slop' in file):
                df_pos_aug = get_df_m6A(path, file)

    for file in directory:
        if 'fasta' in file:
            if ('positive' in file) and ('out' not in file):
                positive_m6A = np.vstack((positive_m6A, m6A_onehot(path, file, df_pos)))
            if ('positive' in file) and ('out' in file):
                positive_m6A = np.vstack((positive_m6A, m6A_onehot(path, file, df_pos_aug)))
    positive_m6A = np.delete(positive_m6A, 0, 0)
    return positive_m6A


def getonehotencoded_files_neg(path, folder):
    directory = os.listdir(path)
    df = pd.DataFrame()
    neg_m6A = np.zeros((1, 400, 5))

    for file in directory:
        if ((file.__contains__('.miclip.filt.bed.out')) & (file.__contains__(folder))):
            df = get_df_m6A(path, file)
    for file in directory:
        if ((file.__contains__(folder)) & (file.__contains__('fa'))):
            neg_m6A = m6A_onehot(path, file, df)
    return neg_m6A


def prepare_raw_dataset_m6A(path, current_folder, augmentation: bool):
    print('Starting preparation of the dataset for %s' % current_folder)
    folder = current_folder
    directory2 = os.listdir(f'{path}{folder}')

    seq_list_pos = np.zeros((1, 400, 5))
    seq_list_neg1 = np.zeros((1, 400, 5))
    seq_list_neg2 = np.zeros((1, 400, 5))

    for folder2 in directory2:
        if folder2 != "fold-4":  # to be split into training and validation
            seq_list_pos_file, seq_list_neg1_file, seq_list_neg2_file = getonehotencoded_files(
                f'{path}{folder}/{folder2}', augmentation)
            # merge the files belonging to the same label
            seq_list_pos = np.vstack((seq_list_pos, seq_list_pos_file))
            seq_list_neg1 = np.vstack((seq_list_neg1, seq_list_neg1_file))
            seq_list_neg2 = np.vstack((seq_list_neg2, seq_list_neg2_file))

        if folder2 == "fold-4":  # for testing
            seq_list_pos_test, seq_list_neg1_test, seq_list_neg2_test = getonehotencoded_files(
                f'{path}{folder}/{folder2}', augmentation)

    seq_list_pos = np.delete(seq_list_pos, 0, 0)
    seq_list_neg1 = np.delete(seq_list_neg1, 0, 0)
    seq_list_neg2 = np.delete(seq_list_neg2, 0, 0)

    return seq_list_pos, seq_list_neg1, seq_list_neg2, seq_list_pos_test, seq_list_neg1_test, seq_list_neg2_test


def prepare_raw_dataset_m6A_new_neg(path, current_folder, path_negatives):
    print('Starting preparation of the dataset for %s' % current_folder)
    folder = current_folder
    directory2 = os.listdir(path + folder)
    seq_list_pos = np.zeros((1, 400, 5))
    seq_list_neg1 = np.zeros((1, 400, 5))

    for folder2 in directory2:
        path_neg = path_negatives
        if folder2 != "fold-4":  # to be split into training and validation 
            seq_list_pos_file = getonehotencoded_files_pos(f'{path}/{folder}/{folder2}')
            seq_list_neg1_file = getonehotencoded_files_neg(path_neg+folder, folder2)

            seq_list_neg1 = np.vstack((seq_list_neg1, seq_list_neg1_file))
            seq_list_pos = np.vstack((seq_list_pos, seq_list_pos_file))

        if folder2 == "fold-4":  # for testing
            seq_list_pos_test = getonehotencoded_files_pos(f'{path}/{folder}/{folder2}')
            seq_list_neg1_test = getonehotencoded_files_neg(path_neg+folder, folder2)
            
    seq_list_pos = np.delete(seq_list_pos, 0, 0)
    seq_list_neg1 = np.delete(seq_list_neg1, 0, 0)

    return seq_list_pos, seq_list_neg1, seq_list_pos_test, seq_list_neg1_test


def prepare_original_dataset(path, current_folder):
    folder = current_folder
    directory2 = os.listdir(f'{path}{folder}')
    seq_list_positive = np.zeros((1, 400, 4))
    seq_list_negative = np.zeros((1, 400, 4))
    seq_list_positive_test = np.zeros((1, 400, 4))
    seq_list_negative_test = np.zeros((1, 400, 4))

    for folder2 in directory2:
        directory3 = os.listdir(f'{path}{folder}/{folder2}')

        if folder2 != "fold-4":  # to be split into training and validation
            seq_list_positive_file = np.zeros((1, 400, 4))
            seq_list_negative_file = np.zeros((1, 400, 4))

            for file in directory3:

                if "fasta" in file and 'out' not in file:
                    fasta_sequences = SeqIO.parse(open(f'{path}{folder}/{folder2}/{file}'),
                                                  'fasta')
                    seq_array = fastatoarray(fasta_sequences)

                    if "positive" in file:
                        seq_list_positive_file = np.vstack((seq_list_positive_file,
                                                            seq_array))

                    if "negative-1" in file:  
                        seq_list_negative_file = np.vstack((seq_list_negative_file,
                                                            seq_array))

            seq_list_negative_file = np.delete(seq_list_negative_file, 0, 0)
            seq_list_positive_file = np.delete(seq_list_positive_file, 0, 0)

            seq_list_positive = np.vstack((seq_list_positive,
                                           seq_list_positive_file))
            seq_list_negative = np.vstack((seq_list_negative,
                                           seq_list_negative_file))

        if folder2 == "fold-4":  # this folder will be used for testing

            seq_list_positive_test_file = np.zeros((1, 400, 4))
            seq_list_negative_test_file = np.zeros((1, 400, 4))

            for file in directory3:

                if "fasta" in file and 'out' not in file:
                    fasta_sequences = SeqIO.parse(open(f'{path}{folder}/{folder2}/{file}'),
                                                  'fasta')
                    seq_array = fastatoarray(fasta_sequences)

                    if "positive" in file:
                        seq_list_positive_test_file = np.vstack((seq_list_positive_test_file,
                                                                 seq_array))

                    if "negative-1" in file:
                        seq_list_negative_test_file = np.vstack((seq_list_negative_test_file,
                                                                 seq_array))

            seq_list_negative_test_file = np.delete(seq_list_negative_test_file, 0, 0)
            seq_list_positive_test_file = np.delete(seq_list_positive_test_file, 0, 0)

            seq_list_positive_test = np.vstack((seq_list_positive_test,
                                                seq_list_positive_test_file))
            seq_list_negative_test = np.vstack((seq_list_negative_test,
                                                seq_list_negative_test_file))

            seq_list_negative_test = np.delete(seq_list_negative_test, 0, 0)
            seq_list_positive_test = np.delete(seq_list_positive_test, 0, 0)

    seq_list_negative = np.delete(seq_list_negative, 0, 0)
    seq_list_positive = np.delete(seq_list_positive, 0, 0)

    # Preparation of the labels
    labels_positive = np.ones((np.shape(seq_list_positive)[0], 1))
    labels_negative = np.zeros((np.shape(seq_list_negative)[0], 1))
    labels_positive_test = np.ones((np.shape(seq_list_positive_test)[0], 1))
    labels_negative_test = np.zeros((np.shape(seq_list_negative_test)[0], 1))
    print('Shape of labels: \n-positive : ', np.shape(labels_negative),
          '\n-negative : ', np.shape(labels_negative),
          '\n-positive validation : ', np.shape(labels_positive_test),
          '\n-negative validation : ', np.shape(labels_negative_test))

    # Merging datasets
    x = np.vstack((seq_list_positive, seq_list_negative))
    x_test = np.vstack((seq_list_positive_test, seq_list_negative_test))
    y = np.vstack((labels_positive, labels_negative))
    y_test = np.vstack((labels_positive_test, labels_negative_test))

    # Splitting dataset  : test, train and validation sets
    test_size = 0.2
    x_train, x_val, y_train, y_val = train_test_split(x, y,
                                                      test_size=test_size,
                                                      shuffle=True)
    x_test, y_test = shuffle(x_test, y_test, random_state=0)

    y_train = y_train.astype("float32")  # actually useful?
    y_val = y_val.astype("float32")
    y_test = y_test.astype("float32")

    x_train = x_train.astype("float32")
    x_val = x_val.astype("float32")
    x_test = x_test.astype("float32")

    print('Shape of datasets: \n-training set : ', np.shape(x_train),
          '\n-validation set : ', np.shape(x_val),
          '\n-testing set : ', np.shape(x_test))

    # Checking class ditribution in the whole dataset and training set

    print('Label frequencies among the dataset')
    plt.hist(y)
    plt.xticks(range(2))
    plt.title('Label Frequency')
    plt.show()

    plt.hist(y_train)
    plt.xticks(range(2))
    plt.title('Label Frequency training set')
    plt.show()

    plt.hist(y_val)
    plt.xticks(range(2))
    plt.title('Label Frequency validation set')
    plt.show()
    return x_train, y_train, x_val, y_val, x_test, y_test


def filter_m6A_from_array(tot, shape):
    tot_filt = np.zeros((1, 400, shape))
    for array in tot:
        if any(array[:, 4]):
            array = np.expand_dims(array, axis=0)
            if shape == 4:
                array = np.delete(array, 4, 2)
            tot_filt = np.vstack((tot_filt, array))
    tot_filt = np.delete(tot_filt, 0, 0)
    return tot_filt


def settingA(pos, neg, pos_t, neg_t, shape, downsampling: bool):
    '''
        downsampling : [True, False]
    '''
    # positive-a: bound and containing at least an m6a site
    # negative-1-a: un-bound and containing at least an m6a site
    # negative-2-a: bound by other proteins and containing at least an m6a site

    # positive, now filtering just for sequences that contain m6A sites

    pos_filt = filter_m6A_from_array(pos, shape)
    pos_t_filt = filter_m6A_from_array(pos_t, shape)
    neg_filt = filter_m6A_from_array(neg, shape)
    neg_t_filt = filter_m6A_from_array(neg_t, shape)

    # downsampling the negative set to match the positive one
    if downsampling is True:
        if len(pos_filt[:][:][:]) > len(neg_filt[:][:][:]):
            pos_filt = pos_filt[:(len(neg_filt[:][:][:]))]
        else:
            neg_filt = neg_filt[:(len(pos_filt[:][:][:]))]
        if len(pos_t_filt[:][:][:]) > len(neg_t_filt[:][:][:]):
            pos_t_filt = pos_t_filt[:(len(neg_t_filt[:][:][:]))]
        else:
            neg_t_filt = neg_t_filt[:(len(pos_t_filt[:][:][:]))]

    # after the filtering, produce the labels
    labels_pos = np.ones((np.shape(pos_filt)[0], 1))
    labels_pos_t = np.ones((np.shape(pos_t_filt)[0], 1))
    labels_neg = np.zeros((np.shape(neg_filt)[0], 1))
    labels_neg_t = np.zeros((np.shape(neg_t_filt)[0], 1))

    return pos_filt, neg_filt, pos_t_filt, neg_t_filt, labels_pos, labels_neg, labels_pos_t, labels_neg_t 


def settingB(pos, neg, pos_t, neg_t, shape, downsampling: bool):
    '''
        downsampling : [True, False]
    '''
    # positive-a: bound and containing at least an m6a site
    # negative-1-a: un-bound and containing at least an m6a site
    # negative-2-a: bound by other proteins and containing at least an m6a site

    # positive, now filtering just for sequences that contain m6A sites

    pos_filt = filter_one_to_one_m6A_from_array(pos, shape)
    pos_t_filt = filter_one_to_one_m6A_from_array(pos_t, shape)
    neg_filt = filter_one_to_one_m6A_from_array(neg, shape)
    neg_t_filt = filter_one_to_one_m6A_from_array(neg_t, shape)

    # downsampling the negative set to match the positive one
    if downsampling is True:
        if len(pos_filt[:][:][:]) > len(neg_filt[:][:][:]):
            pos_filt = pos_filt[:(len(neg_filt[:][:][:]))]
        else:
            neg_filt = neg_filt[:(len(pos_filt[:][:][:]))]
        if len(pos_t_filt[:][:][:]) > len(neg_t_filt[:][:][:]):
            pos_t_filt = pos_t_filt[:(len(neg_t_filt[:][:][:]))]
        else:
            neg_t_filt = neg_t_filt[:(len(pos_t_filt[:][:][:]))]

    # after the filtering, produce the labels
    labels_pos = np.ones((np.shape(pos_filt)[0], 1))
    labels_pos_t = np.ones((np.shape(pos_t_filt)[0], 1))
    labels_neg = np.zeros((np.shape(neg_filt)[0], 1))
    labels_neg_t = np.zeros((np.shape(neg_t_filt)[0], 1))

    return pos_filt, neg_filt, pos_t_filt, neg_t_filt, labels_pos, labels_neg, labels_pos_t, labels_neg_t 


def filter_one_to_one_m6A_from_array(tot, shape):
    tot_filt = np.zeros((1, 400, shape))
    tot_filt_m6a = np.zeros((1, 400, shape))
    tot_filt_no_m6a = np.zeros((1, 400, shape))
    for array in tot:
        if any(array[:, 4]):
            array = np.expand_dims(array, axis=0)
            if shape == 4:
                array = np.delete(array, 4, 2)
            tot_filt_m6a = np.vstack((tot_filt_m6a, array))
        else:
            array = np.expand_dims(array, axis=0)
            if shape == 4:
                array = np.delete(array, 4, 2)
            tot_filt_no_m6a = np.vstack((tot_filt_no_m6a, array))
    tot_filt_m6a = np.delete(tot_filt_m6a, 0, 0)
    tot_filt_no_m6a = np.delete(tot_filt_no_m6a, 0, 0)

    # 1:1 ratio
    if len(tot_filt_m6a) > len(tot_filt_no_m6a):
        tot_filt_m6a = tot_filt_m6a[:(len(tot_filt_no_m6a))]
    if len(tot_filt_m6a) < len(tot_filt_no_m6a):
        tot_filt_no_m6a = tot_filt_no_m6a[:(len(tot_filt_m6a))]

    tot_filt = np.vstack((tot_filt_no_m6a, tot_filt_m6a))
    tot_filt = np.delete(tot_filt, 0, 0)
    return tot_filt


def finalize_dataset(current_folder, pos, neg, pos_t, neg_t, labels_pos, labels_neg, labels_pos_t, labels_neg_t, path_figure, output_folder):

    print('Shape of labels: \n-positive : ', np.shape(labels_pos),
          '\n-negative : ', np.shape(labels_neg), '\n-positive test : ',
          np.shape(labels_pos_t), '\n-negative test : ', np.shape(labels_neg_t))
    # Merging datasets
    x = np.vstack((pos, neg))
    x_test = np.vstack((pos_t, neg_t))
    y = np.vstack((labels_pos, labels_neg))
    y_test = np.vstack((labels_pos_t, labels_neg_t))

    # Creating the validation set
    test_size = 0.2

    # test set

    x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=test_size, shuffle= True)
    x_test, y_test = shuffle(x_test, y_test, random_state=0)

    y_train = y_train.astype("float32") # actually useful?
    y_val = y_val.astype("float32")
    y_test = y_test.astype("float32")

    x_train = x_train.astype("float32")
    x_val = x_val.astype("float32")
    x_test = x_test.astype("float32")

    print('Shape of datasets: \n-training set : ', np.shape(x_train),
          '\n-validation set : ', np.shape(x_val),
          '\n-testing set : ', np.shape(x_test))

    # checking class ditribution in the whole dataset and training set
    (plot_label_freq(y, y_train, y_val, current_folder)).savefig(f'{path_figure}label_freq_{output_folder}.pdf', format='pdf')
    return x_train, y_train, x_val, y_val, x_test, y_test


# tuning of the model architecture
def build_model(hp):

    METRICS = [
      keras.metrics.TruePositives(name='tp'),
      keras.metrics.FalsePositives(name='fp'),
      keras.metrics.TrueNegatives(name='tn'),
      keras.metrics.FalseNegatives(name='fn'),
      keras.metrics.BinaryAccuracy(name='accuracy'),
      keras.metrics.Precision(name='precision'),
      keras.metrics.Recall(name='recall'),
      keras.metrics.AUC(curve='ROC', name='auroc'),
      keras.metrics.AUC(curve='PR', name='auprc')
    ]

    input_shape = (400, 4)
    # Hyperparameter search
    hp_filters = hp.Int('filters', min_value=10, max_value=60, step=10)
    hp_kernel_size = hp.Int('kernel_size', min_value=10, max_value=60, step=5)
    hp_pool_size = hp.Int('pool_size', min_value=1, max_value=10, step=1)
    hp_strides = hp.Int('strides', min_value=1, max_value=10, step=1)
    hp_learning_rate = hp.Float("lr", min_value=1e-4, max_value=1e-2,
                                sampling="log")
    # hp_padding = hp.Choice('padding', ['valid','same'])
    hp_padding = 'same'
    hp_kernel_initializer = hp.Choice('initializer',
                                      ['random_normal', 'random_uniform'])
    hp_activation = hp.Choice("activation", ["relu", "tanh"])

    # Create model
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape=input_shape))
    model.add(keras.layers.Conv1D(filters=hp_filters, kernel_size=hp_kernel_size,
              kernel_initializer=hp_kernel_initializer,
              activation=hp_activation, input_shape=input_shape))
    model.add(keras.layers.MaxPooling1D(pool_size=hp_pool_size,
                                        strides=hp_strides, padding=hp_padding))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Conv1D(filters=hp_filters,
                                  kernel_size=hp_kernel_size,
                                  kernel_initializer=hp_kernel_initializer,
                                  activation=hp_activation,
                                  input_shape=input_shape))
    model.add(keras.layers.MaxPooling1D(pool_size=hp_pool_size,
                                        strides=hp_strides, padding=hp_padding))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Flatten())  # without parameters
    for i in range(hp.Int("num_layers", 1, 3)):
        model.add(
            layers.Dense(
                # Tune number of units separately.
                units=hp.Int(f"units_{i}", min_value=32, max_value=1024, step=32),
                activation=hp_activation
                )
        )
    if hp.Boolean("dropout"):
        model.add(layers.Dropout(rate=0.25))
    model.add(layers.Dense(10, activation=hp_activation))
    model.add(keras.layers.Dropout(0.6))
    model.add(keras.layers.Dense(1, activation='sigmoid'))

    # Compile model
    model.compile(loss='binary_crossentropy',
                  optimizer=tf.keras.optimizers.Adam(hp_learning_rate),
                  metrics=METRICS)

    return model


def build_model_5(hp):

    METRICS = [
      keras.metrics.TruePositives(name='tp'),
      keras.metrics.FalsePositives(name='fp'),
      keras.metrics.TrueNegatives(name='tn'),
      keras.metrics.FalseNegatives(name='fn'), 
      keras.metrics.BinaryAccuracy(name='accuracy'),
      keras.metrics.Precision(name='precision'),
      keras.metrics.Recall(name='recall'),
      keras.metrics.AUC(curve='ROC', name='auroc'),
      keras.metrics.AUC(curve='PR', name='auprc')
    ]

    input_shape = (400, 5)
    # Hyperparameter search
    hp_filters = hp.Int('filters', min_value=10, max_value=60, step=10)
    hp_kernel_size = hp.Int('kernel_size', min_value=10, max_value=60, step=5)
    hp_pool_size = hp.Int('pool_size', min_value=1, max_value=10, step=1)
    hp_strides = hp.Int('strides', min_value=1, max_value=10, step=1)
    hp_learning_rate = hp.Float("lr", min_value=1e-4, max_value=1e-2,
                                sampling="log")
    # hp_padding = hp.Choice('padding', ['valid','same'])
    hp_padding = 'same'
    hp_kernel_initializer = hp.Choice('initializer',
                                      ['random_normal', 'random_uniform'])
    hp_activation = hp.Choice("activation", ["relu", "tanh"])

    # Create model
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape=input_shape))
    model.add(keras.layers.Conv1D(filters=hp_filters, kernel_size=hp_kernel_size,
              kernel_initializer=hp_kernel_initializer,
              activation=hp_activation, input_shape=input_shape))
    model.add(keras.layers.MaxPooling1D(pool_size=hp_pool_size,
                                        strides=hp_strides, padding=hp_padding))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Conv1D(filters=hp_filters,
                                  kernel_size=hp_kernel_size,
                                  kernel_initializer=hp_kernel_initializer,
                                  activation=hp_activation,
                                  input_shape=input_shape))
    model.add(keras.layers.MaxPooling1D(pool_size=hp_pool_size,
                                        strides=hp_strides, padding=hp_padding))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Flatten())  # without parameters
    for i in range(hp.Int("num_layers", 1, 3)):
        model.add(
            layers.Dense(
                # Tune number of units separately.
                units=hp.Int(f"units_{i}", min_value=32, max_value=1024, step=32),
                activation=hp_activation
                )
        )
    if hp.Boolean("dropout"):
        model.add(layers.Dropout(rate=0.25))
    model.add(layers.Dense(10, activation=hp_activation))
    model.add(keras.layers.Dropout(0.6))
    model.add(keras.layers.Dense(1, activation='sigmoid'))

    # Compile model
    model.compile(loss='binary_crossentropy',
                  optimizer=tf.keras.optimizers.Adam(hp_learning_rate),
                  metrics=METRICS)

    return model


def create_baseline_standard():  # this is the original model architecture --> constant 

    METRICS = [
        keras.metrics.TruePositives(name='tp'),
        keras.metrics.FalsePositives(name='fp'),
        keras.metrics.TrueNegatives(name='tn'),
        keras.metrics.FalseNegatives(name='fn'), 
        keras.metrics.BinaryAccuracy(name='accuracy'),
        keras.metrics.Precision(name='precision'),
        keras.metrics.Recall(name='recall'),
        keras.metrics.AUC(curve='ROC', name='auroc'),
        keras.metrics.AUC(curve='PR', name='auprc')
    ]

    # Parameters
    input_shape = (400, 4)
    filters = 30
    kernel_size = 25
    pool_size = 2
    strides = 2

    # Create model
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape=input_shape))
    model.add(keras.layers.Conv1D(filters=filters, kernel_size=kernel_size,
              kernel_initializer='random_normal',
              activation='relu'))
    model.add(keras.layers.MaxPooling1D(pool_size=pool_size, strides=strides,
                                        padding='valid'))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Conv1D(filters=filters, kernel_size=kernel_size,
              kernel_initializer='random_normal',
              activation='relu', input_shape=input_shape ))
    model.add(keras.layers.MaxPooling1D(pool_size=pool_size, strides=strides,
                                        padding='valid'))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(1024, activation='relu'))
    model.add(keras.layers.Dropout(0.6))
    model.add(keras.layers.Dense(128, activation='relu'))
    model.add(keras.layers.Dropout(0.6))
    model.add(keras.layers.Dense(1, activation='sigmoid'))

    # Compile model
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=METRICS)

    return model


def create_baseline_model(input_shape):  # GOING TO DELETE THIS IF THE DOUBLE SAVING OF THE MODEL WORKS OUT

    METRICS = [
        keras.metrics.TruePositives(name='tp'),
        keras.metrics.FalsePositives(name='fp'),
        keras.metrics.TrueNegatives(name='tn'),
        keras.metrics.FalseNegatives(name='fn'), 
        keras.metrics.BinaryAccuracy(name='accuracy'),
        keras.metrics.Precision(name='precision'),
        keras.metrics.Recall(name='recall'),
        keras.metrics.AUC(curve='ROC', name='auroc'),
        keras.metrics.AUC(curve='PR', name='auprc')
    ]
    # inserting the best model hyperparameters keeping the input shape variable

    # Parameters
    filters = 50
    kernel_size = 10
    pool_size = 10
    strides = 8
    learning_rate = 0.0007744959

    # create model
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape=input_shape))
    model.add(keras.layers.Conv1D(filters=filters, kernel_size=kernel_size,
              kernel_initializer='random_uniform',
              activation='tanh'))
    model.add(keras.layers.MaxPooling1D(pool_size=pool_size, strides=strides, padding='same'))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Conv1D(filters=filters, kernel_size=kernel_size,
              kernel_initializer='random_uniform',
              activation='tanh'))
    model.add(keras.layers.MaxPooling1D(pool_size=pool_size, strides=strides, padding='same'))
    model.add(keras.layers.Dropout(0.3))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(768, activation='tanh', use_bias=True, kernel_initializer='glorot_uniform'))
    model.add(keras.layers.Dense(256, activation='tanh', use_bias=True, kernel_initializer = 'glorot_uniform'))
    model.add(keras.layers.Dense(288, activation='tanh', use_bias=True,kernel_initializer='glorot_uniform'))
    model.add(keras.layers.Dense(10, activation='tanh', use_bias=True, kernel_initializer='glorot_uniform'))
    model.add(keras.layers.Dropout(0.6))
    model.add(keras.layers.Dense(1, activation='sigmoid', use_bias=True, kernel_initializer='glorot_uniform'))

    # Compile model
    model.compile(loss='binary_crossentropy', optimizer=tf.keras.optimizers.Adam(learning_rate), metrics=METRICS)
    return model


def training_model(model, class_weights, x_train, y_train, x_val, y_val, x_test, y_test, folder, resultspath, output_folder, path_figure):
    # should include the name of the folder in the path ( for example best_models_baseline)
    '''
    output_folder is the folder but not the actual name of the file: that will be given idividually to every model depending on the RBP in exam--> there is one output_folder for loop of same type of model. . Example : _model_settingA or _model_baseline
    List of output_folder avilable : 
    '''

    my_callbacks = [
        tf.keras.callbacks.EarlyStopping(patience=5), #better to reduce to 2 for the search 
        #tf.keras.callbacks.ModelCheckpoint(filepath=f'{resultspath}/model_checkpoints/model.{val_loss:.2f}.h5'),
        #tf.keras.callbacks.TensorBoard(log_dir='./logs'),
        ]

    # fitting the model
    history = model.fit(x_train, y_train, epochs=200, validation_data=(x_val,y_val), class_weight=class_weights, callbacks=my_callbacks)

    # Generate figure
    (plot_metrics_training(history, folder)).savefig(f'{path_figure}metrics_training{output_folder}.pdf', format='pdf')

    model.save(f'{resultspath}best_models/{output_folder}/%s' % folder)


def testing_model(path, current_folder, x_test, y_test, predictionspath, output_folder, path_figure):

    reconstructed_model = keras.models.load_model(path)
    y_pred = reconstructed_model.predict(x_test)
    threshold = 0.5
    y_pred_thresh = np.where(y_pred > threshold, 1, 0)

    # saving the true and predicted labels
    df = pd.DataFrame(y_pred, columns=['Predictions'])
    df['True'] = y_test

    # Generate figure
    (plot_metrics_prediction(y_test, y_pred, y_pred_thresh, current_folder)).savefig(f'{path_figure}metrics_predictions_{output_folder}.pdf', format='pdf')

    # saving the predictions (0,1) in a csv file -- in probabilities
    df.to_csv(f"{predictionspath}{output_folder}/%s.csv" % current_folder)
    print("Prediction saved.")


def prepare_best_model(x_train, y_train, x_val, y_val, x_test, y_test, resultspath, current_folder):

    model = build_model(keras_tuner.HyperParameters())

    tuner = keras_tuner.RandomSearch(
        hypermodel=build_model,
        objective="val_accuracy",
        max_trials=200,
        executions_per_trial=1,
        overwrite=True,  # 'False' because I found the best combination
        directory="my_dir",
        project_name="models",
        max_consecutive_failed_trials=10,
    )

    tuner.search(x_train, y_train, epochs=5, validation_data=(x_val, y_val))

    # makes sense to save two models, one for 4 channels one for 5 channels
    # Get the top 2 models
    #models = tuner.get_best_models(num_models=4)
    #best_model = models[0]
    # Build the model.
    #best_model.build(input_shape=(None, 400, 5))
    #best_model.summary() #save to be used later 

    best_hps = tuner.get_best_hyperparameters(5)
    print(best_hps[0])
    # SAVE THESE HYPERPARAMETERS AND USE THEM FOR create_baseline_model()
    #at this point the optimal parameters have been identified, the next part of training and prediction is the same for all models 

    # Build the model with the best hp for both input sizes --> will it work?
    #models = tuner.get_best_models(num_models=4)
    #best_model = models[0]
    model = build_model(best_hps[0]) # just if the input shape is not inside the fuction 
    model_5 = build_model_5(best_hps[0])
    #model = best_model.build(input_shape=(None, 400, 4))
    #model_5 = best_model_5.build(input_shape=(None, 400, 5))

    # saving the best model for testing and aggregated evaluation later
    model.save(f'{resultspath}best_models/best_hps_model{current_folder}')
    model_5.save(f'{resultspath}best_models/best_hps_model_5{current_folder}')

def plot_metrics_prediction(y_test, y_pred, y_pred_thresh, folder):
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 30))
    fig.suptitle('Metrics Training among the dataset %s' % folder)
    false_pr, true_pr, thresholds = metrics.roc_curve(y_test, y_pred, drop_intermediate=False)

    ax[0].plot(false_pr, true_pr, label=" (AUC=%.2f)" % metrics.roc_auc_score(y_test, y_pred))
    ax[0].plot([0, 1], [0, 1], color="grey", label="Random Classifier", linestyle="--")
    ax[0].set_xlabel("False Positive Rate")
    ax[0].set_ylabel("True Positive Rate")
    ax[0].set_ylim(0, 1)
    ax[0].set_xlim(0, 1)
    ax[0].grid(color="#CCCCCC")
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].legend()

    precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
    from sklearn.metrics import average_precision_score

    
    #plot = metrics.PrecisionRecallDisplay.from_predictions(y_test, y_pred)
    #ax[1].plot(plot)
    ax[1].plot(recall, precision, color='purple')
    ax[1].set_title('Precision-Recall Curve')
    ax[1].set_ylabel('Precision')
    ax[1].set_xlabel('Recall')
    ax[1].legend(['AC=%.2f' % metrics.average_precision_score(y_test, y_pred)])
    
    matrix = confusion_matrix(y_test, y_pred_thresh)
    names = ['True Neg', 'False Pos', 'False Neg', 'True Pos']
    counts = ['{0:0.0f}'.format(value) for value in matrix.flatten()]
    percentages = ['{0:.2%}'.format(value) for value in matrix.flatten()/np.sum(matrix)]
    labels = [f'{v1}\n{v2}\n{v3}' for v1, v2, v3 in zip(names, counts, percentages)]
    labels = np.asarray(labels).reshape(2, 2)
    ax[2] = sns.heatmap(matrix, annot=labels, fmt='', cmap='Blues')

    return fig


def plot_metrics_training(history, folder):
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 30))

    # summarize history for accuracy
    axs[0].set_title('model accuracy %s' % folder)
    axs[0].plot(history.history['accuracy'])
    axs[0].plot(history.history['val_accuracy'])
    axs[0].set_ylabel('accuracy')
    axs[0].set_xlabel('epoch')
    axs[0].legend(['train', 'test'], loc='upper left')

    # summarize history for loss
    axs[1].set_title('model loss %s' % folder)
    axs[1].plot(history.history['loss'])
    axs[1].plot(history.history['val_loss'])
    axs[1].set_ylabel('loss')
    axs[1].set_xlabel('epoch')
    axs[1].legend(['train', 'test'], loc='upper left')

    # summarize history for auroc
    axs[2].set_title('model auroc and auprc %s' % folder)
    axs[2].plot(history.history['auroc'])
    axs[2].plot(history.history['val_auroc'])
    axs[2].plot(history.history['auprc'])
    axs[2].plot(history.history['val_auprc'])
    axs[2].set_ylabel('auroc and auprc')
    axs[2].set_xlabel('epoch')
    axs[2].legend(['auroc_train', 'auroc_test', 'auprc_train', 'auprc_test'], loc='upper left')

    fig.suptitle('Best model evaluation for %s' % folder)
    return fig


def plot_label_freq(y, y_train, y_val, current_folder):
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 30))
    fig.suptitle('Label frequencies among the dataset %s' % current_folder)
    ax[0].hist(y)
    ax[0].set_xticks(range(2))
    ax[0].set_title('Label Frequency')

    ax[1].hist(y_train)
    ax[1].set_xticks(range(2))
    ax[1].set_title('Label Frequency training set')

    ax[2].hist(y_val)
    ax[2].set_xticks(range(2))
    ax[2].set_title('Label Frequency validation set')
    return fig

def create_table_latex(data_table):
    data_table[['RBP', 'Cell line']] = data_table["rbp"].apply(lambda x: pd.Series(str(x).split("_")))
    data_table = data_table.drop(columns = ['rbp'])

    cols = data_table.columns.tolist()
    cols = ['Cell line', 'RBP', 'roc_auc_score', 'pr_auc_score']
    data_table = data_table[cols]

    print((((data_table.sort_values('RBP')).groupby('Cell line', group_keys=False).apply(lambda x: x))).to_latex(index=False,
                      formatters={"name": str.upper},
                      float_format="{:.4f}".format,
    )) 
    
def boxplot_roc(df, title):
    plt.figure(figsize=(6,6))
    sns.boxplot(x=df["roc_auc_score"], color = '#48AAAD' ).set(title= title)
    plt.xlabel('AUROC scores')
    return plt

def boxplot_pr(df, title):
    plt.figure(figsize=(6,6))
    sns.boxplot(x=df["pr_auc_score"], color = '#F9812A' ).set(title= title)
    plt.xlabel('PR scores')
    return plt

def boxplot_double_roc(df1, df2, title):
    plt.figure(figsize=(6,6))
    data_auc = pd.concat([df1, df2])
    sns.boxplot(x=data_auc["roc_auc_score"], y = data_auc['model'], hue = data_auc['model']).set(title= title)
    #sns.boxplot(x=df["roc_auc_score"], color = '#48AAAD' ).set(title= title)
    plt.xlabel('AUROC scores')
    return plt


def boxplot_double_pr(df1, df2, title):
    plt.figure(figsize=(6,6))
    data_auc = pd.concat([df1, df2])
    sns.boxplot(x=data_auc["pr_auc_score"], y = data_auc['model'], hue = data_auc['model']).set(title= title)
    #sns.boxplot(x=df["roc_auc_score"], color = '#48AAAD' ).set(title= title)
    plt.xlabel('PR scores')
    return plt