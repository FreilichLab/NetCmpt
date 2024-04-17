# -*- coding: utf-8 -*-
import pandas as pd
import pickle
import sys
import os
import networkx as nx
import numpy as np

#Input are

#tab delmited lists separted by space

#genomes.tsv: tab delmited lists of ECs separted by space

def get_essential_compounds_number(ind,df_ann, env, df_compounds,compound_header):

    compound_query_set =  set(df_ann.loc[ind,env]) 

    compound_ref_set = set(df_compounds[compound_header].tolist())

    set_inter = compound_query_set.intersection(compound_ref_set)

    essential_list = list(set_inter)

    return len(essential_list)

def get_compounds_number(ind,df_ann, env):
    return(len(df_ann.loc[ind,env]))

def get_set_diffrence(ind,df_ann, env,df_gen_env):

    #get index x

    ind_x = df_ann.loc[ind,'index_x']

    #get index y

    ind_y = df_ann.loc[ind,'index_y']

    #get index x boolan (row) in genomes_env  dataframe

    ind_x_bool = df_gen_env['index'] == ind_x

    #get index y boolan (row) in genomes_env  dataframe

    ind_y_bool = df_gen_env['index'] == ind_y
    
    #get get seeds in row x as dataframe

    set_x= set(df_gen_env[ind_x_bool][env].tolist()[0])

    #get get seeds in row y as dataframe

    set_y = set(df_gen_env[ind_y_bool][env].tolist()[0])

    #get  set diffrence set_x - set_y

    set_x_y = set_x.difference(set_y)

    return(list(set_x_y))

def extract_seeds(EClist):
        try:
             with open(BaseFolder+"database/DB.pickle", 'rb') as handle:
                 DB = pickle.load(handle)
        except:
            DB = pd.read_pickle(BaseFolder+"data/DB/DB.pickle")      
          
        df_el = DB['full_enzymes_labels_jun.txt']
        df_ecMapping = DB['ec_reac_mapping_jun.txt']
        df_reactions = DB['reactions_3_balanced.txt']
        df_ec_to_compoundIndex = DB['compound_labels_jun.txt']
        rlist = df_ecMapping["Reactions"].values.tolist()   
        df_ecMapping["Reactions"] = [i.split(":")[1].lstrip().rstrip().split(" ") for i in rlist]    
        IndexEnzymeList = df_el[["Index", "Enzyme_Code"]].values.tolist()
        DictEL = {}
        
        for i in IndexEnzymeList:
            DictEL[i[1]] = i[0]
        
        #fix list strings to list
        #edgeR_row = ast.literal_eval(edgeR_row)
        
        #EClist = edgeR_row[0]
        print("Extracting EClist")
    
        EClist_ = []
        
        for i in EClist:
            try:
                EClist_.append(DictEL[i])
            except:
                _=""
                
        df_ecMapping = df_ecMapping[df_ecMapping["Index"].isin(EClist_)]
        ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
        flat_list = []
                
        for i in ListOfRelevantReactions:
            for j in i:
                try:
                    flat_list.append(int(j))
                except:
                    _=""
                    
        df_reactions = df_reactions[df_reactions['Reaction_index'].isin(flat_list)]
        #l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Left"].values.tolist()]
        #r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Reaction_Right"].values.tolist()]
        l = df_reactions["Reaction_Left"].values.tolist()
        r = df_reactions["Reaction_Right"].values.tolist()     
        DictComp = {}
        
        for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
            DictComp[i[0]] = i[1]        
        
        df_tmp = pd.DataFrame()
        df_tmp["l"]=l
        df_tmp["r"]=r
        df_tmp = df_tmp.explode("l").explode("r")
        df_tmp_grouped = pd.DataFrame({'count' : df_tmp.groupby( ["l","r"] ).size()}).reset_index()
    
        G_list = [ (i[0], i[1], i[2]) for i in df_tmp_grouped.values.tolist() ]
        G_seeds = nx.DiGraph()
        G_seeds.add_weighted_edges_from(G_list) 
        seeds = [len(c) for c in sorted(nx.strongly_connected_components_recursive(G_seeds), key=len)]    
        seeds += list(set(list(np.setdiff1d(df_tmp["l"].values,df_tmp["r"].values))))    
        MySeeds = [str(DictComp[int(i)]) for i in seeds]
        CompList = list(set(df_tmp["l"].values.tolist()+df_tmp["r"].values.tolist()))
        CompList = [str(DictComp[int(i)]) for i in CompList]
        MySeeds = list(set(MySeeds))
        CompList = list(set(CompList))
        EClist = list(set(EClist))
        print("Seeds were calculated")
        return MySeeds, CompList, EClist

def is_consist(LeftSides, RightSides, directions, BucketOfCompounds, Reaction_Direction):
    IsConsist = []
    new = BucketOfCompounds
    for row in range(len(LeftSides)):

        if ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '>'))) or \
                ((set(LeftSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):

            IsConsist.append(1)
            new = new + RightSides[row]
            # print(str(LeftSides[row])+" IS IN "+str(new))

        elif ((set(RightSides[row]).issubset(set(new))) and ((Reaction_Direction[row] == '='))):

            IsConsist.append(1)
            new = new + LeftSides[row]
            # print(str(LeftSides[row])+" IS IN "+str(new))

        else:
            IsConsist.append(0)
            # print(str(LeftSides[row])+" is NOT in "+str(new))

    return IsConsist, new

def simulation(df_genomes_envs_, resfolder, analysisMode_, DB,genomes_,envs_):
    # Open log file
    file_object = open(resfolder + 'Simulation_log.txt', 'a')

    #procees raw input from extranal files
    #organisms = input1[0].values.tolist()
    #organisms = [i.split(" ")[0] for i in organisms]
    #enzymes = input1[0].values.tolist()
    #enzymes = [i.split(" ")[1:] for i in enzymes]

    organisms = df_genomes_envs_['index'].tolist()

    enzymes = df_genomes_envs_[genomes_].tolist()

    # Create all pairs##############################
    # Present multiple analysis or a single analysis for each env
    if (analysisMode_ == "interaction"):
        orgs = []
        enzs = []

        for i, vali in enumerate(organisms):
            for j, valj in enumerate(organisms):
                if (i < len(organisms)) and (j > i):
                    # print(vali+" "+valj)
                    orgs.append(str(vali) + "-" + str(valj))
                    enzs.append(list(set(enzymes[i] + enzymes[j])))

        organisms = organisms + orgs
        enzymes = enzymes + enzs

    file_object.write('#### List of organisms #### \n')
    file_object.write(str(organisms) + ' \n')
    file_object.write('#### List of enzymes #### \n')
    file_object.write(str(enzymes) + ' \n')

    input1_ = pd.DataFrame()
    input1_["Organism"] = organisms
    input1_["Enzymes_List"] = enzymes

    # envs = list(input2["Env_Content"].values)
    # envs = [i.strip().split(" ")[1:] for i in envs]

    #rawinput

    #envs = input2[0].values.tolist()


    #envs = [i.strip().split(" ")[1:] for i in envs]

    envs = df_genomes_envs_[envs_].tolist()

    file_object.write('#### List of Env #### \n')
    file_object.write(str(envs) + ' \n')

    ## PATCH - FIX INPUT2##
    D = {}
    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        D[str(i[1]).strip()] = str(i[0]).strip()

    for i, vali in enumerate(envs):
        for j, valj in enumerate(vali):
            envs[i][j] = D[str(valj).strip()]

    #####

    # Dictionary - relevant enzymes of each organism
    OrgEnzDict = {}

    # Dictionary - Compound index to Compound Code
    CompoundDict = {}

    for i in DB['compound_labels_jun.txt'][["Index", "Compound_Code"]].values.tolist():
        CompoundDict[i[0]] = i[1]

    # enviornments:
    listOfEnvsOfOrganisms = []
    newOrganismsCompounds = []

    SimSteps_df = pd.DataFrame(columns=["Cycle", "Organism", "InitialEnv", "numberOfCompounds", "Compounds"])
    simulation_steps = []

    for count_, k in enumerate(envs):
        file_object.write('#### Simulation - env ' + str(k) + ' \n')
        file_object.write(str(k) + ' \n')

        # For each Organism:
        Enzymes_List = []
        k = [int(kk) for kk in k]
        newOrganismsCompounds = []
        for i, vali in enumerate(input1_["Organism"].values):
            Enzymes_List = input1_[input1_["Organism"] == vali]
            Enzymes_List = list(Enzymes_List["Enzymes_List"].values)[0]
            Enzymes_List = [j for j in Enzymes_List if str(j) != 'nan']

            file_object.write('#### Simulation - Organism ' + str(vali) + ' \n')
            file_object.write(str(vali) + ' \n')

            # Enzymes_list code to index:
            ll = DB["full_enzymes_labels_jun.txt"][["Index", "Enzyme_Code"]].values.tolist()
            d = {}
            for j in ll:
                d[j[1]] = j[0]

            Enzymes_List_ = []
            for j in Enzymes_List:
                try:
                    Enzymes_List_.append(d[j])
                except:
                    # print(str(j)+" Not in list")
                    file_object.write('#### Simulation - Enzyme list : ' + str(j) + ' Not in list \n')

            # Find reactions consist input Compounds: check if all compounds in envirnment (k) are in the same side of the equation
            # AND an enzyme (j) is in the list --> TAKE the other side of equation according to equilibrium > < =

            df = DB["ec_reac_mapping_jun.txt"][DB["ec_reac_mapping_jun.txt"]['Index'].isin(Enzymes_List_)]

            RelevantReactions = df["Reactions"].values.tolist()
            RelevantReactions = [x.lstrip().rstrip().split(":")[1].lstrip().rstrip().split(" ") for x in
                                 RelevantReactions]
            RelevantReactions = [item for sublist in RelevantReactions for item in sublist]

            for x, valx in enumerate(RelevantReactions):
                try:
                    RelevantReactions[x] = int(valx)
                except:
                    # print("cant make "+valx+" int")
                    file_object.write('#### Simulation - Relevant Reactions : ' + 'cant make ' + str(valx) + ' int  \n')

            RelevantReactions = [x for x in RelevantReactions if x != '']

            OrgEnzDict[i] = RelevantReactions

            df = DB["reactions_3_balanced.txt"][
                DB["reactions_3_balanced.txt"]['Reaction_index'].isin(RelevantReactions)]

            RL = df["Reaction_Left"].values.tolist()
            RR = df["Reaction_Right"].values.tolist()
            Reaction_Direction_ = df["Reaction_Direction"].values.tolist()

            OnlyNewCompounds = []
            newCompounds = []
            prevCompounds = k
            C = 0

            IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds,
                                          Reaction_Direction=Reaction_Direction_)
            simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])

            OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)))

            while set(newCompounds) != set(prevCompounds):
                OnlyNewCompounds = list(set(newCompounds).intersection(set(prevCompounds)))
                simulation_steps.append([count_, C, vali, k, len(set(newCompounds)), list(set(newCompounds))])
                prevCompounds = list(set(newCompounds))
                IC, newCompounds = is_consist(RL, RR, Reaction_Direction_, prevCompounds,
                                              Reaction_Direction=Reaction_Direction_)
                C += 1

            file_object.write('#### Simulation : ' + str(C) + ' Cycles \n')

            # Here build final network snapshot (input cytoscape)
            # 1 - Filter by relevant enzyme (done) 2 - filter by relevant reactions (only reactions which are subsets of the final env)
            # 3 - flatten the list

            newCompoundsCodes = [CompoundDict[k] for k in newCompounds]

            if C > 0:
                newOrganismsCompounds.append(list(set(newCompoundsCodes)))
            else:
                newOrganismsCompounds.append([])

        listOfEnvsOfOrganisms.append(newOrganismsCompounds)

    SimSteps_df = pd.DataFrame(simulation_steps,
                               columns=["EnvIndex", "Cycle", "Organism", "InitialEnv", "numberOfCompounds",
                                        "Compounds"])
    Final_df = pd.DataFrame(listOfEnvsOfOrganisms, columns=organisms)

    # Close the log file
    file_object.close()
    return Final_df, SimSteps_df


try:
    BaseFolder = sys.argv[1]
    organismsFileName = sys.argv[2]
    analysisMode = sys.argv[3]
    output_file = sys.argv[4]
    UserInputFiles = [organismsFileName]
    # print(sys.argv)
    # print(UserInputFiles)
except:
    print("Error in executing. Using default parameters")
    BaseFolder = "./"
    UserInputFiles = [BaseFolder + "genomes.tsv"]
    analysisMode = "Single-species"
    output_file = BaseFolder + 'compete.tsv'
try:
    with open(BaseFolder + "database/DB.pickle", 'rb') as handle:
         DB = pickle.load(handle)
except:
    DB = pd.read_pickle(BaseFolder + "database/DB.pickle")

 # 1 - Take Input2 compounds (medium)
Input1 = pd.read_csv(UserInputFiles[0], header=None, sep = '\t', names = ['species', 'ECs'])

Input1.set_index('species', drop = False, inplace = True)

#convert space delimieted list to a python list

#df_EC_list

df_EC_list = Input1['species'].apply(lambda x: Input1.loc[x, 'ECs'].split(' ')).to_frame()

df_EC_list.columns = ['ECs']

df_EC_list.reset_index(inplace = True)

df_EC_list.set_index('species', drop = False, inplace = True)

#df_seeds

df_extract_seeds = df_EC_list['species'].apply(lambda x: extract_seeds(df_EC_list.loc[x,'ECs'])).to_frame()

# Add index and spcecies  columns to EClist

df_EC_list.drop(columns = 'species', inplace = True)

df_EC_list.reset_index(inplace = True)

df_EC_list.reset_index(inplace = True)

#Extract seed using extract_seeds from netcom

#df_seeds
#
df_extract_seeds.columns = ['extract_seeds']
#
df_extract_seeds.reset_index(inplace = True)
#
df_extract_seeds.set_index('species', drop = False, inplace = True)
#
##put each result of extract seeds in separte columns (
## split list of lists into multiple columns
df_extract_seeds = pd.DataFrame(df_extract_seeds['extract_seeds'].tolist(), index = df_extract_seeds.index)
#
df_extract_seeds.columns = ['seeds', 'CompList', 'EClist']
#
df_seeds = df_extract_seeds.drop(columns = ['CompList', 'EClist'])

df_seeds.reset_index(inplace = True)
#
#join genomes and seeed (envs) to a single df

#From this calculate matablic capabilities in in optimal enviroment (seeds)

df_genomes_envs = pd.merge(left = df_EC_list, right = df_seeds, how = 'inner', on = ['species'])

##create cross_prodcut of X, Y _indexes 
#
#df1 = pd.DataFrame({'X':[1,2]})
#df2 = pd.DataFrame({'Y':[1,2]})
#df1_2_cross = df1.merge(df2, how='cross')

df_genomes_cross = df_genomes_envs[['index']].merge(df_genomes_envs[['index']], how = 'cross')

df_species_cross = df_genomes_envs[['species']].merge(df_genomes_envs[['species']], how = 'cross')

df_species_cross.reset_index(inplace = True)

df_gen_env_cross = df_genomes_cross.merge(df_genomes_envs,
        how = 'left', left_on = 'index_x', right_on = 'index')

df_gen_env_cross.drop(columns = ['index'], inplace = True)

df_gen_env_cross.reset_index(inplace = True)

df_test_cross_diff = df_gen_env_cross['index'].apply(get_set_diffrence,df_ann = df_gen_env_cross,env = 'seeds',df_gen_env= df_genomes_envs).to_frame()             

df_test_cross_diff.columns = ['envs_diffrence']

df_test_cross_diff.reset_index(inplace = True)

df_gen_env_cross = df_gen_env_cross.merge(df_test_cross_diff, how = 'inner', on = ['index'] )
##get essential compounds
df_essential_compounds = DB['color_index.txt']['C_Code'].to_frame()

# get num_essential_seeds_compounds

df_test_comp = df_genomes_envs['index'].apply(get_essential_compounds_number,df_ann = df_genomes_envs,env = 'seeds',df_compounds = df_essential_compounds,compound_header = 'C_Code').to_frame()             

df_test_comp.columns = ['num_essential_seeds_compounds']

df_test_comp.reset_index(inplace = True)

df_genomes_envs = df_genomes_envs.merge(df_test_comp, on = ['index'], how = 'inner')

# get num_seeds_compounds

df_test_comp = df_genomes_envs['index'].apply(get_compounds_number,df_ann = df_genomes_envs,env = 'seeds').to_frame()             

df_test_comp.columns = ['num_seeds_compounds']

df_test_comp.reset_index(inplace = True)

df_genomes_envs = df_genomes_envs.merge(df_test_comp, on = ['index'], how = 'inner')

def run_simulation_final(ind,df_ann,genomes, envs):

    df_final, df_simsteps = simulation(df_genomes_envs_ = df_ann.iloc[[ind],:],resfolder=BaseFolder,
            analysisMode_=analysisMode, DB=DB, genomes_ = genomes,envs_ = envs)

    return df_final.loc[0,ind]

df_temp = df_genomes_envs.copy()['index'].apply(run_simulation_final,df_ann = df_genomes_envs.copy(deep = True),
        genomes = 'ECs', envs = 'seeds').to_frame() 


df_temp.columns = ['metabolic_capacity_compounds']

df_temp.reset_index(inplace = True)

# Add 'metabolic_capacity_compounds' genomes_envs

df_genomes_envs = df_genomes_envs.merge(df_temp, on = ['index'], how = 'inner')

#### compute number of essential metabolic_capacity_compounds

df_test_comp = df_genomes_envs['index'].apply(get_essential_compounds_number,
df_ann = df_genomes_envs,env = 'metabolic_capacity_compounds',df_compounds = df_essential_compounds, compound_header = 'C_Code').to_frame()             

df_test_comp.columns = ['num_essential_metabolic_capacity_compounds']

df_test_comp.reset_index(inplace = True)

df_genomes_envs = df_genomes_envs.merge(df_test_comp, on = ['index'], how = 'inner')

#### compute number of metabolic_capacity_compounds

df_test_comp = df_genomes_envs['index'].apply(get_compounds_number,df_ann = df_genomes_envs,env = 'metabolic_capacity_compounds').to_frame()             

df_test_comp.columns = ['num_metabolic_capacity_compounds']

df_test_comp.reset_index(inplace = True)

df_genomes_envs = df_genomes_envs.merge(df_test_comp, on = ['index'], how = 'inner')

#### run simulaton on cross (pairwise)

df_temp = df_gen_env_cross['index'].apply(run_simulation_final,df_ann = df_gen_env_cross,
        genomes = 'ECs', envs = 'envs_diffrence').to_frame() 

df_temp.columns = ['metabolic_capacity_compounds']

df_temp.reset_index(inplace = True)

# Add 'metabolic_capacity_compounds' genomes_envs

df_gen_env_cross = df_gen_env_cross.merge(df_temp, on = ['index'], how = 'inner')

#### compute number of essential metabolic_capacity_compounds

df_test_comp = df_gen_env_cross['index'].apply(get_essential_compounds_number,
df_ann = df_gen_env_cross,env = 'metabolic_capacity_compounds',df_compounds = df_essential_compounds, compound_header = 'C_Code').to_frame()

df_test_comp.columns = ['num_essential_metabolic_capacity_compounds']

df_test_comp.reset_index(inplace = True)

df_gen_env_cross = df_gen_env_cross.merge(df_test_comp, on = ['index'], how = 'inner')

#### compute number of metabolic_capacity_compounds

df_test_comp = df_gen_env_cross['index'].apply(get_compounds_number,df_ann = df_gen_env_cross,env = 'metabolic_capacity_compounds').to_frame()

df_test_comp.columns = ['num_metabolic_capacity_compounds']

df_test_comp.reset_index(inplace = True)

df_gen_env_cross = df_gen_env_cross.merge(df_test_comp, on = ['index'], how = 'inner')

def compute_EMO(ind,df_gen_env_inter,df_gen_env_cross_inter,df_emo_inter,df_mo_inter):

    #get index x
    print(ind)
    ind_x = df_gen_env_cross_inter.loc[ind,'index_x']


    #get index y

    ind_y = df_gen_env_cross_inter.loc[ind,'index_y']

    #get index x boolan (row) in df_gen_env_cross_inter

    ind_x_bool_cross = df_gen_env_cross_inter['index'] == ind_x

    #get index y boolan (row) in df_gen_env_cross_inter

    ind_y_bool_cross = df_gen_env_cross_inter['index'] == ind_y

    #get index x boolan (row) in df_gen_env_inter

    ind_x_bool = df_gen_env_inter['index'] == ind_x

    #get index y boolan (row) in df_gen_env_cross_inter

    ind_y_bool = df_gen_env_inter['index'] == ind_y

    #find number of compounds (fs)

    f_emo = df_gen_env_inter[ind_x_bool]['num_essential_metabolic_capacity_compounds']

    f_mo = df_gen_env_inter[ind_x_bool]['num_metabolic_capacity_compounds']

    f_emo_diff = df_gen_env_cross_inter.loc[ind,'num_essential_metabolic_capacity_compounds']

    f_mo_diff = df_gen_env_cross_inter.loc[ind,'num_metabolic_capacity_compounds']

    if int(f_emo) == 0:
        emo = 0
        print("f_emo = 0 for ind_x: ", ind_x, 'ind_y: ', ind_y)
    else:
        emo = 1 - int(f_emo_diff)/int(f_emo)

    if int(f_mo) == 0:
        mo = 0
        print("f_mo = 0 for ind_x: ", ind_x, 'ind_y: ', ind_y)
    else:

        mo = 1 - int(f_mo_diff)/int(f_mo)

    df_emo_inter.loc[ind_x, ind_y] = emo

    df_mo_inter.loc[ind_x, ind_y] = mo

    return [emo,mo]

num_species = len(df_genomes_envs['index']) 

#Effictive metabolic overlap (EMO)

df_emo = pd.DataFrame(columns = list(range(num_species)), index = list(range(num_species)))

#Metabolic Overlap (MO)

df_mo = pd.DataFrame(columns = list(range(num_species)), index = list(range(num_species)))

df_emo_mo = df_gen_env_cross['index'].apply(compute_EMO,df_gen_env_inter=df_genomes_envs, df_gen_env_cross_inter = df_gen_env_cross,
        df_emo_inter = df_emo,df_mo_inter = df_mo).to_frame()

df_emo_mo.columns = ['emo_mo']
       
df_emo_mo[['emo','mo']] = pd.DataFrame(df_emo_mo['emo_mo'].tolist(), index = df_emo_mo.index)

df_emo_mo.drop(columns = ['emo_mo'], inplace = True)

df_emo_mo.reset_index(inplace = True)

df_species_cross = df_species_cross.merge(df_emo_mo,how = 'inner', on =['index'])

df_species_cross.to_csv(output_file, sep = '\t', index = False)
