#-------------------------------------------------------------------------------
# Name:        IMPACCT 
# Purpose:     Calculate quality index scores 
#       
# Author:      Sean P. Sheehy
#
# Created:     18/05/2016; Version 2.5
# Copyright:   (c) Sean P. Sheehy 2016
# Licence:     
#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
sns.set(font_scale=1.8)

# Function to calculate Hellinger distance metric from mean and variance of 
# control and treatment groups.  Returns a distance value in the range of [0-1]
# overlap parameter indicates whether degree of overlap between distributions
# is returned, or degree of non-overlap (default)
def Hellinger_Uni(mu1,sigma1,mu2,sigma2, overlap=False):
    if sigma1 == 0.0 or sigma2 == 0.0:
        return 0.0
    else:
        left = np.sqrt((2.0*sigma1*sigma2)/(sigma1**2 + sigma2**2))
        right = np.exp(-0.25*((mu1 - mu2)**2/(sigma1**2 + sigma2**2)))
        
        if overlap:
            return (left*right)
        else:
            return (1.0 - (left*right))


# Function to calculate the Strictly Standardized Mean Difference from mean and
# variance of control and treatment groups used to determine direction of 
# therapeutic effect
def SSMD(mu1,sigma1,mu2,sigma2):
    return ((mu2 - mu1)/np.sqrt((sigma1**2 + sigma2**2)))


# Function to calculate IMPACCT score for an individual parameter using the
# mean and standard deviations of the control condition and the condition under
# comparison against the control condition. scale parameter allows scoring range
# to be converted from [0-1] to desired range through multiplication 
def Parameter_Score(con_mu, con_sigma, exp_mu, exp_sigma, effect=False, overlap=False, scale=100):
    # If direction of therapeutic effect is provided, then use Strictly Standardized
    # Mean Difference (SSMD) to determine direction of difference and calculate 
    # IMPACCT score as Hellinger distance value with sign of SSMD, otherwise
    # output absolute value of Hellinger distance scores
    if effect:
        # Calculate the Hellinger distance between the control and experimental 
        # condition using mean and standard deviation values
        hellinger = Hellinger_Uni(con_mu, con_sigma, exp_mu, exp_sigma, overlap)
    
        # Calculate the Strictly Standardized Mean Difference (SSMD) between control
        # and experimental condition to determine the direction of the difference
        # between the two populations using mean and standard deviation values
        ssmd_sign = np.sign(SSMD(con_mu, con_sigma, exp_mu, exp_sigma))
        
        score = hellinger*ssmd_sign
        
    else:
        # Calculate the Hellinger distance between the control and experimental 
        # condition using mean and standard deviation values
        score = Hellinger_Uni(con_mu, con_sigma, exp_mu, exp_sigma, overlap)
        
    # Return the IMPACCT score scaled according to scale input parameter
    return score*scale


# Main function of IMPACCT scoring routine that takes as input a file containing
# either columns of raw data or mean and standard deviation values for parameters
# to be compared, a file containing the 'expected therapeutic direction', and
# the name of the condition that is to be used as the 'gold standard' or control
# condition that all other conditions are to be compared against
def IMPACCT(data_file, raw_data=False, effect=False, effect_file=None, overlap=False, out_rem=False, standard='Control', scale=100):
    
    ########## READ INPUT DATA FROM FILES #########
    # Read data from input file into DataFrames (input files can be Excel
    # workbooks with data organized into sheets, or .csv flat files)
    # Data should be organized such that columns represent parameters
    # and rows represent measurements from individual subjects
    if data_file.lower().endswith('.csv'):
        input_data = pd.read_csv(data_file)
    elif data_file.lower().endswith('.xlsx'):
        input_data = pd.read_excel(data_file, sheetname=None)
    else:
        print data_file, "is in an unrecognized file format"    
    
    # Read expected therapeutic effect mappings for each parameter used for
    # IMPACCT score from Excel worksheets into a DataFrame
    # 0 = no expected change, 1 = increase expected, -1 = decrease expected
    if effect:
        if effect_file.lower().endswith('.csv'):
            expected_effect = pd.read_csv(effect_file)
        elif data_file.lower().endswith('.xlsx'):
            expected_effect = pd.read_excel(effect_file, sheetname='Sheet1')
        else:
            print effect_file, "is in an unrecognized file format"
    
    ########## CALCULATE IMPACCT SCORES FROM INPUT DATA #########
    # Empty DataFrame to hold IMPACCT scores calculated for each parameter
    # in each experimental condition
    IMPACCT_scores = pd.DataFrame(index=input_data.keys(), columns=input_data[standard].keys())
    
    # If input file contains raw data values, then calculate mean and standard
    # deviation before computing parameter scores, otherwise just compute scores
    # from mean and std dev values read in from file for each parameter
    if raw_data:
        # If removal of outliers is desired, remove values from each parameter
        # column that are greater than 3 standard deviations from the mean
        if out_rem:
            input_data = input_data[(np.abs(stats.zscore(input_data)) < 3)]
            
        # Calculate mean and standard deviation values for each parameter in the
        # 'gold standard' or control condition to compare against
        control = {}
        for parameter in input_data[standard].keys():
            control[parameter] = [input_data[standard][parameter].mean(), input_data[standard][parameter].std()]
    
        # Calculate mean and standard deviation values for each parameter in the
        # conditions that are to be tested against the control condition
        for condition in input_data.keys():
            if condition != standard:
                for parameter in input_data[condition].keys():
                    avg = input_data[condition][parameter].mean()
                    std = input_data[condition][parameter].std()
                    
                    # Calculate IMPACCT score for each parameter, comparing
                    # against the value in the control condition
                    score = Parameter_Score(control[parameter][0], control[parameter][1], avg, std, effect=effect, overlap=overlap, scale=scale)
                    # Check to make sure the sign of the score matches the
                    # expected therapeutic effect, and adjust sign of output
                    # score if necessary
                    if effect:
                        score = score*expected_effect[parameter].values[0]
                    # Save score to IMPACCT_scores DataFrame for graphing
                    IMPACCT_scores.set_value(condition, parameter, score)
    else:
        # Calculate mean and standard deviation values for each parameter in the
        # 'gold standard' or control condition to compare against
        control = {}
        for parameter in input_data[standard].keys():
            control[parameter] = [input_data[standard][parameter].values[0], input_data[standard][parameter].values[1]]
    
        # Calculate mean and standard deviation values for each parameter in the
        # conditions that are to be tested against the control condition
        for condition in input_data.keys():
            if condition != standard:
                for parameter in input_data[condition].keys():
                    avg = input_data[condition][parameter].values[0]
                    std = input_data[condition][parameter].values[1]
                    
                    # Calculate IMPACCT score for each parameter, comparing
                    # against the value in the control condition
                    score = Parameter_Score(control[parameter][0], control[parameter][1], avg, std, effect=effect, overlap=overlap, scale=scale)
                    # Check to make sure the sign of the score matches the
                    # expected therapeutic effect, and adjust sign of output
                    # score if necessary
                    if effect:
                        score = score*expected_effect[parameter].values[0]
                    # Save score to IMPACCT_scores DataFrame for graphing
                    IMPACCT_scores.set_value(condition, parameter, score)
        
    
    # Remove the empty first row of the IMPACCT_scores DataFrame corresponding
    # to the control condition and replace NaNs with 0 values
    #IMPACCT_scores = IMPACCT_scores.ix[1:].fillna(0)
    IMPACCT_scores = IMPACCT_scores.ix[1:]
    
    # Calculate combined IMPACCT score from all parameters for each condition
    # and add as a new column to the IMPACCT_scores DataFrame
    IMPACCT_scores['Combined'] = IMPACCT_scores.mean(axis=1)
    
    # Replace NaN values from data set so that the heatmap function won't complain
    # and flag them for masking
    IMPACCT_scores = IMPACCT_scores.fillna(-1)
    
    #print IMPACCT_scores.mean(axis=1)
    
    # Generate heatmap to illustrate IMPACCT score outcome
    # Generate a custom diverging colormap
    if effect:
        cmap = sns.diverging_palette(10, 220, center='dark', as_cmap=True)
        vmin=scale*-1
    else:
        cmap='RdYlGn'
        vmin=0
        
    # Draw a heatmap with the numeric values in each cell
    sns.heatmap(IMPACCT_scores, mask=IMPACCT_scores < 0, vmin=vmin, vmax=scale, cmap=cmap, annot=True, fmt=".1f", linewidths=.5)
    plt.rc('font',weight='bold', size=32)
    sns.plt.show()
    
    # Return DataFrame with IMPACCT scores 
    return IMPACCT_scores
                     

def main():
    pass

if __name__ == '__main__':

    # Input data set in Excel file with tabs indicating each condition
    input_data = 'IMPACCT_Test_Input.xlsx'
    
    # Expected therapeutic effect for wound healing parameters
    #wound_effect = 'IMPACCT_Wound_Effect.xlsx'

    # Calculate IMPACCT scores for all of the experimental conditions and parameters
    IMPACCT_Output = IMPACCT(input_data, raw_data=False, standard='Control', effect=False, overlap=True, out_rem=False, scale=100)
    # Save IMPACCT scores to file
    IMPACCT_Output.to_csv('IMPACCT_Test_Output.csv')
    
    
    # Test case for Hellinger distance function to ensure correct computation
    sanity = Hellinger_Uni(139.0,48.0,146.0,51.0)