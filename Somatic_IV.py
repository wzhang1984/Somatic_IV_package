# Version 1.1: added the option to include covariates in the 1st 3 steps: 1. network-regularized Lasso, 2. OLS, 3. Aalen hazard model

import glob
import warnings
import os
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from scipy.stats import t
from sklearn.linear_model import ElasticNetCV, LassoCV, Lasso, LassoLarsIC
import statsmodels.stats.multitest as multi
import statsmodels.api as sm
import statsmodels.formula.api as smf
import multiprocessing
from io import StringIO

import networkx as nx

class MR_Egger_model():
    def __init__(self, Dir, disease, genes_in_network, IVs_in_network, PPR, gene2P0, X, G0, G0_prime, metadata, 
                 delay_entry=False, LASSO_method='CV', covariates=[]):
        self.Dir = Dir
        self.disease = disease
        self.genes_in_network = genes_in_network
        self.IVs_in_network = IVs_in_network
        self.PPR = PPR
        self.gene2P0 = gene2P0
        self.X = X
        self.G0 = G0
        self.G0_prime = G0_prime
        self.metadata = metadata
        self.delay_entry = delay_entry
        self.LASSO_method = LASSO_method
        self.covariates = covariates
        
    # This function calculates each driver's weights based on 
    # its network propagation score from the target gene 
    def get_driver_weights(self):
        P0 = self.gene2P0.loc[self.genes_in_network, :]

        P = np.dot(P0, self.PPR)
        P = pd.DataFrame(P, index=P0.index, columns=P0.columns)

        W = P.loc[:, self.IVs_in_network]

        # Normalize by the highest number of each row
        W = (W.T/W.apply(lambda row: row.nlargest(1).values[-1], axis=1)).clip(upper=1).T
        
        # # Only keep interactions within 2 hops
        # mask_2hops = W.copy()*0
        # for gene in mask_2hops.index:
        #     mask_2hops.loc[gene, mask_2hops.columns.isin(nx.single_source_shortest_path_length(self.Graph, gene, cutoff=2).keys())]=1
        # W = np.multiply(W, mask_2hops)

        W = W.dropna()
        return W
    
    # This function calculates the linear models Pr~IV
    def run_OLS(self, x, G0):
        lm_stats = []
        for feature in G0.columns:
            G = pd.concat([G0.loc[:, [feature]], self.metadata.loc[:, self.covariates]], axis=1, join='inner')
            G['intersect'] = 1
            x = x.loc[G.index]
            md = sm.OLS(x, G)
            mdf = md.fit()
            lm_stats.append(pd.read_html(StringIO(mdf.summary().tables[1].as_html()), header=0, index_col=0)[0].iloc[[0]])
            lm_stats[-1]['F-statistic'] = pd.read_html(StringIO(mdf.summary().tables[0].as_html()))[0].iloc[2, 3]
        lm_stats = pd.concat(lm_stats)
        return lm_stats

    # This function performs MR-Egger regression based on summary statistics
    def MR_Egger_regression(self, lm_XG, lm_YG, Regulator):
        BetaYG = lm_YG.loc[:, 'coef']
        BetaXG = lm_XG.loc[:, 'coef']
        seBetaYG = lm_YG.loc[:, 'std err']
        seBetaXG = lm_XG.loc[:, 'std err']
        
        BYG = BetaYG*np.sign(BetaXG) 
        BXG = abs(BetaXG)
        BXG = sm.add_constant(BXG)
        
        md = sm.WLS(BYG, BXG, weights=1/seBetaYG**2)
        mdf = md.fit()
        
        MREggerBeta0 = mdf.params.iloc[0]
        MREggerBeta1 = mdf.params.iloc[1]
        SE0 = mdf.bse.iloc[0]/min(1, np.sqrt(mdf.scale))
        SE1 = mdf.bse.iloc[1]/min(1, np.sqrt(mdf.scale))
        DF = len(BetaYG)-2
        
        MRBeta0_p = 2*(1-t.cdf(abs(MREggerBeta0/SE0),DF))
        MRBeta1_p = 2*(1-t.cdf(abs(MREggerBeta1/SE1),DF))
        MRBeta0_CI = MREggerBeta0 + np.array([-1,1])*t.ppf(0.975, df=DF)*SE0
        MRBeta1_CI = MREggerBeta1 + np.array([-1,1])*t.ppf(0.975, df=DF)*SE1

        MREggerResults = pd.DataFrame([MREggerBeta0, MREggerBeta1, SE0, SE1, MRBeta0_CI[0], MRBeta0_CI[1], MRBeta1_CI[0], MRBeta1_CI[1], 
                                       MREggerBeta0/SE0, MREggerBeta1/SE1, MRBeta0_p, MRBeta1_p]).T
        MREggerResults.columns = ['Beta_0E', 'Beta_E', 'Beta_0E_SE', 'Beta_E_SE', 'Beta_0E_025', 'Beta_0E_975', 'Beta_E_025', 'Beta_E_975', 
                                             'Beta_0E_t', 'Beta_E_t', 'Beta_0E_p', 'Beta_E_p']
        MREggerResults['Regulator'] = Regulator
        MREggerResults['nIVs'] = len(BXG)
        MREggerResults['IVs'] = ', '.join([i for i in (BXG['coef'].sort_values()[::-1]).index])
        MREggerResults['IVs_coefs'] = ', '.join(['{:.2f}'.format(i) for i in (BXG['coef'].sort_values()[::-1]).values])
        
        return MREggerResults
    
    # Run MR-Egger regression for a single protein across all phenotypes
    def run_MR_Egger(self, Regulator):
        MREggerResults = pd.DataFrame()

        # 2. BetaXG: Regulator-IV Associations (OLS)
        lm_XG = self.lm_XG.loc[self.lm_XG['Regulator']==Regulator]

        # 3. BetaYG: Survival-IF Associations (Aalen)
        lm_YG = self.lm_YG

        # 4. MR-Egger regression
        index = lm_XG.index.intersection(lm_YG.index)
        if len(index) >= 3:
            MREggerResults = self.MR_Egger_regression(lm_XG.loc[index], lm_YG.loc[index], Regulator)

        return MREggerResults
    
    # 1. IV selection: Pr-IV Associations (multivariate model)
    # This function calculates the LASSO models Protein_Activity ~ alterations
    def run_LASSO(self):
        warnings.simplefilter(action='ignore')

        preds = []
        models = []
        scores = []
        for gene in self.W.index:
            w = self.W.loc[gene]
            G =  self.G0 * w
            G_prime = self.G0_prime * w

            G = pd.concat([G, self.metadata.loc[:, self.covariates]], axis=1, join='inner')
            G_prime = pd.concat([G_prime, self.metadata.loc[:, self.covariates]], axis=1, join='inner')
            if len(self.covariates) != 0:
                G_prime.loc[:, self.covariates] = 0
            x = self.X.loc[G.index, gene]

            if self.LASSO_method == 'CV':
                regr = LassoCV(cv=5, random_state=0, n_jobs=5) # CV
            elif self.LASSO_method == 'BIC':
                regr = LassoLarsIC(criterion='bic') # BIC
            regr.fit(G, x)

            # Check F-statistic 
            score = regr.score(G_prime.loc[x.index], x)
            n = len(x)
            if len(self.covariates) == 0:
                k = sum(regr.coef_!=0)
            else:
                k = sum(regr.coef_[:-len(self.covariates)]!=0)
            F = score/(1-score)*(n-k-1)/k

            models.append(list(regr.coef_))
            models[-1].append(regr.intercept_)
            
            preds.append(regr.predict(G_prime))

            if self.LASSO_method == 'CV':
                MSE = regr.mse_path_.mean(1).min() # CV
                scores.append([gene, MSE, score, F, k]) # CV
            elif self.LASSO_method == 'BIC':
                scores.append([gene, score, F, k]) # BIC

        models = pd.DataFrame(models)
        models.index = self.W.index
        columns = list(G.columns)
        columns.append('intercept')
        models.columns = columns
        
        preds = pd.DataFrame(preds)
        preds.index = self.W.index
        preds.columns = G_prime.index

        scores = pd.DataFrame(scores)
        scores = scores.set_index(0)
        if self.LASSO_method == 'CV':
            scores.columns = ['MSE', 'score', 'F', 'k'] # CV
        elif self.LASSO_method == 'BIC':
            scores.columns = ['score', 'F', 'k'] # BIC    

        return models, preds, scores
    
    # 2. BetaXG: Pr-IV Associations
    def run_lm_XG(self, Regulator):
        x = self.X.loc[:, [Regulator]].dropna()
        # IV selection by LassoCV
        IV2coef_Lasso = self.models.loc[Regulator].drop(['intercept']+self.covariates)
        IVs = IV2coef_Lasso.loc[IV2coef_Lasso!=0].index
        G0 = self.G0.loc[x.index, IVs]
        lm_XG = self.run_OLS(x, G0)
        lm_XG['Regulator'] = Regulator
        return lm_XG
        
    # 3. BetaYG: Survival~IV Associations
    # This function run Aalen additive hazard model
    def run_Aalen(self):
        # run Aalen additive hazard model in R
        self.train.to_csv('tmp/train.csv')
        self.metadata_lite.to_csv('tmp/metadata.csv')
        if self.delay_entry:
            os.system("Rscript run_Aalen_model_delay_entry.R") # TEMPUS OS
        else:
            os.system("Rscript run_Aalen_model.R") # TEMPUS PFI
        lm_YG = pd.read_csv('tmp/results.txt', sep='\t', index_col=0)
        lm_YG = lm_YG.loc[:, ['Coef.', 'Robust SE', 'z', 'P-val']]
        lm_YG.columns = ['coef', 'std err', 'z', 'P-val']
        return lm_YG
            
    def fit(self, n=10, load_lasso_models=False, load_BetaXG=False, load_BetaYG=False, LASSO_only=False):        
        with multiprocessing.Pool(processes=n) as pool:
            # W is a regulator-by-IV matrix. 
            # W_ij represent the random walk score from regulator_i to the IV_j that carries somatic alterations
            print('# Calculating IV weights')
            self.W = self.get_driver_weights()
            self.W.to_csv('data/{}/{}_IV_weights.csv'.format(self.Dir, self.disease))
        
            print('# 1. IV selection: Regulator-IVs (multivariate model)')
            if load_lasso_models:
                self.preds = pd.read_csv('data/{}/{}_LASSO_preds.csv'.format(self.Dir, self.disease), index_col=0)
                self.scores = pd.read_csv('data/{}/{}_LASSO_scores.csv'.format(self.Dir, self.disease), index_col=0)
                self.models = pd.read_csv('data/{}/{}_LASSO_models.csv'.format(self.Dir, self.disease), index_col=0)
            else:
                self.models, self.preds, self.scores = self.run_LASSO()
                self.models.to_csv('data/{}/{}_LASSO_models.csv'.format(self.Dir, self.disease))
                self.preds.to_csv('data/{}/{}_LASSO_preds.csv'.format(self.Dir, self.disease))
                self.scores.to_csv('data/{}/{}_LASSO_scores.csv'.format(self.Dir, self.disease))
            print(self.scores.head())
            if self.LASSO_method == 'CV':
                self.Regulators_explainable = self.scores.loc[(self.scores['MSE']<=0.99) & (self.scores['score']>=0.01) & 
                                                              (self.scores['k']>=3)].index # CV
            elif self.LASSO_method == 'BIC':
                self.Regulators_explainable = self.scores.loc[(self.scores['score']>=0.01) & (self.scores['k']>=3)].index # BIC
            
            if LASSO_only:
                pool.close()
                pool.join()
                return
            
            print('# 2. BetaXG: Regulator-IV Associations (univariate OLS)')
            if load_BetaXG:
                self.lm_XG = pd.read_csv('data/{}/{}_OLS_XG.csv'.format(self.Dir, self.disease), index_col=0)
            else:
                # Only include proteins whose activity can be explained by the genetic alterations
                self.lm_XG = pool.map(self.run_lm_XG, list(self.Regulators_explainable))
                self.lm_XG = pd.concat(self.lm_XG)
            print(self.lm_XG.head())
            self.lm_XG.to_csv('data/{}/{}_OLS_XG.csv'.format(self.Dir, self.disease))
            
            print('# Build input data for the Aalen model')
            
            self.train = self.G0_prime.copy()
            # Option to include factors (e.g. age, sex, race) in the metadata
            if self.delay_entry:
                self.metadata_lite = self.metadata.loc[:, ['event', 'time_to_event', 'time_delay_entry']+self.covariates] # TEMPUS OS
            else:
                self.metadata_lite = self.metadata.loc[:, ['event', 'time_to_event']+self.covariates] # TEMPUS PFI, TCGA
            print('# 3. BetaYG: Survival-IV Associations (univariate Aalen model)')
            if load_BetaYG:
                self.lm_YG = pd.read_csv('data/{}/{}_Aalen_TG.csv'.format(self.Dir, self.disease), index_col=0)
            else:
                self.lm_YG = self.run_Aalen()
            self.lm_YG = self.lm_YG.dropna()
            print(self.lm_YG.head())
            self.lm_YG.to_csv('data/{}/{}_Aalen_TG.csv'.format(self.Dir, self.disease))
            
            print('# 4. MR-Egger regression')
            self.MREggerResults = pool.map(self.run_MR_Egger, list(self.Regulators_explainable))
            self.MREggerResults = pd.concat(self.MREggerResults).sort_values('Beta_E_p')
            self.MREggerResults['Beta_E_fdr'] = multi.multipletests(self.MREggerResults['Beta_E_p'], method='fdr_bh')[1]
            self.MREggerResults = self.MREggerResults.set_index('Regulator')
            self.MREggerResults.to_csv('data/{}/{}_MREggerResults.csv'.format(self.Dir, self.disease))
            print(self.MREggerResults.head())
            pool.close()
            pool.join()
        return self.MREggerResults