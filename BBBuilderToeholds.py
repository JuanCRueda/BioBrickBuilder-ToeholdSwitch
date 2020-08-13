import pandas as pd
from subprocess import Popen, PIPE
from sklearn.linear_model import LogisticRegression
import numpy as np
from random import randint
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import streamlit as st
import base64
from io import BytesIO
from dna_features_viewer import GraphicFeature, GraphicRecord

def main():

        st.cache(persist=True)
        def Promoter_Selection(strength):
        	'''Loads the Anderson promoters' data as a pd.DataFrame, with a given strength computes the difference with the reported strength, checks for compability with the selected standards and returns the Id, Sequence and relative strength as dict'''
        	promoter_df=pd.read_csv('Andersonpromoters.csv')
        	promoter_df=promoter_df.dropna(axis=0,how='any')
        	promoter_df['Distance']=abs(promoter_df['Measured Strengthb']-strength)
        	while True:
        		promoter_index=promoter_df['Distance'].idxmin(axis=1)
        		if test_standard(promoter_df.loc[promoter_index,'Sequencea'],enzyme_set):
        			promoter_data={'id':promoter_df.loc[promoter_index,'Identifier'],'seq':promoter_df.loc[promoter_index,'Sequencea'],'Strength':promoter_df.loc[promoter_index,'Measured Strengthb']}
        			return promoter_data
        		else:
        			promoter_df=promoter_df.drop(promoter_index,axis=0)

        st.cache(persist=True)
        def RBS_Selection(strength):
        	'''Loads a data set with data from selected BioBrick RBS as a pd.DataFrame, with a given strength computes the difference with the reported strength, checks for compability with the selected standards and returns the Id, Sequence and relative strength as dict'''
        	RBS_df=pd.read_csv('RBS.csv')
        	RBS_df=RBS_df.dropna(axis=0,how='any')
        	RBS_df['Distance']=abs(RBS_df['Strength']-strength)
        	while True:
        		RBS_index=RBS_df['Distance'].idxmin(axis=1)
        		if test_standard(RBS_df.loc[RBS_index,'Sequence'],enzyme_set):
        			RBS_data={'id':RBS_df.loc[RBS_index,'Identifier'],'seq':RBS_df.loc[RBS_index,'Sequence'],'Strength':RBS_df.loc[RBS_index,'Strength']}
        			return RBS_data
        		else:
        			RBS_df=RBS_df.drop(RBS_index)

        def MFE_Toehold(seq):
        	'''Return the MFE (kcal/mol) of the RNA secondary strucure of a given sequence'''
        	MFE_calc=Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE)
        	Result=MFE_calc.communicate(seq.encode())
        	return float(Result[0][-9:-3])

        def MFE_Hybridization(seq_toehold,seq_target):
        	'''Return the MFE (kcal/mol) of the hybridization of two given RNA sequences'''
        	MFE_calc=Popen('C:\Program Files (x86)\ViennaRNA Package\RNAcofold.exe', stdin=PIPE, stdout=PIPE)
        	Input='>Seq_toehold\n'+seq_toehold+'\n>Seq_target\n'+seq_target
        	Result=MFE_calc.communicate(Input.encode())
        	return float(Result[0][-9:-3])

        @st.cache(persist=True)
        def train_classifier(treshold):
            toehold_df=pd.read_excel('Toehold_Data_Processed.xlsx')
            toehold_df['Class']=np.where(toehold_df['On/Off ratio']>=treshold, 1,0)
            X_values=toehold_df.loc[:,'MFE Toehold':'MFE Hybridization']
            y_values=toehold_df['Class']
            clf=LogisticRegression().fit(X_values, y_values)
            return clf
        
        #@st.cache(persist=True)
        def ToeholdSequence_gen(seq,pool,rbs,treshold):
        	'''From a given DNA sequences generates a pool of n random subsequences, checks them for standard compability and specificity, assembles the toehold switch according to the specification by Green et al., 2014, calculares the MFE of the secondary structure and hybridation, predicts is it would be over the minimum (treshold) On/Off ratio and returns the best candidate'''
        	seq=seq.lower()
        	targets=[]
        	toeholds=[]
        	reversed_targets=[]
        	mutated_targets=[]
        	for _ in range(pool):
        		n=randint(0,len(seq)-30)
        		target=seq[n:n+30]
        		target_dna=Seq(target,generic_dna)
        		target_reverse=str(target_dna.reverse_complement())
        		target_mutated=target[:6]+'ATG'+target[9:18]
        		target_assembled=target_reverse+'CAAG'+rbs['seq']+target_mutated+'AACCTGGCGGCAGCGCAAAAG'
        		if test_standard(target_assembled,enzyme_set):
        			if BLAST_test(target,organism):
        				targets.append(target)
        				toeholds.append(target_assembled)
        				reversed_targets.append(target_reverse)
        				mutated_targets.append(target_mutated)
        	toehold_df=pd.DataFrame({'Target':targets,'Toehold':toeholds,'Reversed':reversed_targets,'Mutated':mutated_targets})
        	Toehold_values=[]
        	Hibrid_values=[]
        	for ind in list(toehold_df.index.values):
        		Toehold_values.append(MFE_Toehold(toehold_df.loc[ind,'Toehold']))
        		Hibrid_values.append(MFE_Hybridization(toehold_df.loc[ind,'Toehold'],toehold_df.loc[ind,'Target']))
        	toehold_df['MFE Toehold']=Toehold_values
        	toehold_df['MFE Hybridization']=Hibrid_values
        	x=toehold_df.loc[:,'MFE Toehold':'MFE Hybridization']
        	clf=train_classifier(treshold)
        	y_pred=clf.predict_proba(x)
        	probabilities=[c[1] for c in y_pred]
        	toehold_df['class']=probabilities
        	toehold_index=toehold_df['class'].idxmin(axis=1)
        	toehold={'Toehold_seq':toehold_df.loc[toehold_index,'Toehold'],'seq_target':toehold_df.loc[toehold_index,'Target'],'id_reversed':'NA','seq_reversed':toehold_df.loc[toehold_index,'Reversed'],'id_spacer':'NA','seq_spacer':'CAAG','id_mutated':'NA','seq_mutated':toehold_df.loc[toehold_index,'Mutated'],'id_linker':'NA','seq_linker':'AACCTGGCGGCAGCGCAAAAG'}
        	return toehold
        
        @st.cache(persist=True)
        def get_standard(standard_list):
        	'''Create a set containing the prohibited restriction sites (subsequences) from a givin standard list'''
        	enzyme_set=set()
        	if 'RFC10' in standard_list:
        		enzyme_set.update(['gaattc','tctaga','actagt','ctgcag','gcggccgc'])
        	if 'RFC12' in standard_list:
        		enzyme_set.update(['gaattc','actagt','gctagc','ctgcag','gcggccgc','cagctg','ctcgag','tctaga','gctcttc','gaagagc'])
        	if 'RFC21' in standard_list:
        		enzyme_set.update(['gaattc','agatct','ggatcc','ctcgag'])
        	return enzyme_set

        def test_standard(seq,enzyme_set):
        	'''Given a DNA sequence and a set of prohibited restriction site, checks for the prohibited sites within the sequence, returns True if no prohibited sites where detected'''
        	seq=seq.lower()
        	for enz in enzyme_set:
        		if enz in seq:
        			return False
        	return True

        def BLAST_test(seq,organism):
        	'''Runs a BLAST to test wether a given sequence is present at a given organism, returns False if it's present'''
        	organism=organism.lower()
        	result_handle=NCBIWWW.qblast('blastn','nt',seq)
        	blast_record=NCBIXML.read(result_handle)
        	for alignment in blast_record.alignments:
        		for hsp in alignment.hsps:
        			if int(hsp.identities)==len(seq):
        				seq_title=alignment.title
        				seq_title=seq_title.lower()
        				if organism in seq_title:
        					return False
        	return True

        @st.cache(persist=True)
        def Check_DNA(seq):
            '''Check if a given sequence is a valid DNA sequence, returns True if it is a valid DNA sequence'''
            seq=seq.lower()
            for base in seq:
            	if base not in 'atcg':
            		return False
            return True

        @st.cache(persist=True)
        def load_reporters(standard_list):
            ''' Loads a selected reporter proetin dataset as pd.DataFrame, checks or standard compliances and drop the noncompliant instances, returns a pd.DataFrame with the compliant instances'''
            standard_list=list(standard_list)
            reporter_df=pd.read_excel('Reporters.xlsx')
            reporter_df=reporter_df.drop_duplicates('Sequence')
            for ind in list(reporter_df.index.values):
            	Not_complient=False
            	for standard in standard_list:
            		if standard not in reporter_df.loc[ind,'Standard']:
            			Not_complient=True
            			break
            		if Not_complient:
            			reporter_df.drop(ind)
            return reporter_df

        @st.cache(persist=True)
        def Assembler(promoter,rbs,toehold,output_df):
        	'''Assemble the BioBrick from the given parts, returns the assembled sequence and the annotation table'''
        	output_df=output_df.reset_index()
        	types=['BioBrick prefix','Promoter',"Trigger'",'Spacer','RBS','Toehold complement','Linker','Output','Terminator','BioBrick suffix']
        	seq='GAATTCGCGGCCGCTTCTAGAG'+promoter['seq']+'TACTAGAG'+toehold['Toehold_seq']+output_df.loc[0,'Sequence']+'TACTAGAG'+'tcacactggctcaccttcgggtgggcctttctgcgtttatatactagagagagaatataaaaagccagattattaatccggcttttttattattt'+'TACTAGTAGCGGCCGCTGCAG'
        	ids=['BioBrick prefix',promoter['id'],toehold['id_reversed'],toehold['id_spacer'],rbs['id'],toehold['id_mutated'],toehold['id_linker'],output_df.loc[0,'Id'],'BBa_B0014','BioBrick suffix']
        	Starts=[1,24]
        	Ends=[23]
        	MoreInfo=['NA','Strength: '+str(promoter['Strength']),'NA','Green et al., 2014','Strength: '+str(rbs['Strength']),'NA','Green et al., 2014',output_df.loc[0,'Description'],'Double E.coli terminator, widely used','NA']
        	Sequences=['GAATTCGCGGCCGCTTCTAGAG',promoter['seq'],toehold['seq_reversed'],toehold['seq_spacer'],rbs['seq'],toehold['seq_mutated'],toehold['seq_linker'],output_df.loc[0,'Sequence'],'tcacactggctcaccttcgggtgggcctttctgcgtttatatactagagagagaatataaaaagccagattattaatccggcttttttattattt','TACTAGTAGCGGCCGCTGCAG']
        	Ends.append(Starts[-1]+len(promoter['seq']))
        	Starts.append(Ends[-1]+9)
        	Ends.append(Starts[-1]+len(toehold['seq_reversed']))
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+len(toehold['seq_spacer']))
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+len(rbs['seq']))
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+len(toehold['seq_mutated']))
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+len(toehold['seq_linker']))
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+len(output_df.loc[0,'Sequence']))
        	Starts.append(Ends[-1]+9)
        	Ends.append(Starts[-1]+95)
        	Starts.append(Ends[-1]+1)
        	Ends.append(Starts[-1]+21)
        	assembled_df=pd.DataFrame({'Type':types,'BioBrick Id':ids,'Start':Starts,'End':Ends,'More Information':MoreInfo,'Sequence':Sequences})
        	return seq, assembled_df

        @st.cache(persist=True)
        def to_excel(df):
            output = BytesIO()
            writer = pd.ExcelWriter(output, engine='xlsxwriter')
            df.to_excel(writer, sheet_name='Sheet1')
            writer.save()
            processed_data = output.getvalue()
            return processed_data

        @st.cache(persist=True)
        def get_table_download_link(df):
                """Generates a link allowing the data in a given panda dataframe to be downloaded
                in:  dataframe
                out: href string
                """
                val = to_excel(df)
                b64 = base64.b64encode(val)  # val looks like b'...'
                return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="Toeholdswitch.xlsx">Download as Excel</a>'

        st.title('BioBrick Builder: Toeholdswitch designer')
        st.sidebar.title('BioBrick Builder: Toeholdswitch designer')
        st.sidebar.markdown('Fill out the fields with your preferences for your toeholdswitch part')
        st.sidebar.subheader('For which organism is your part?')
        organism=st.sidebar.selectbox('Organism',('Escherichia coli','Saccharomyces cerevisiae','Insects','Arabidopsis thaliana'))
        st.sidebar.markdown('')
        st.sidebar.subheader('Standards')
        standards=st.sidebar.multiselect('Select the standards that you part should be compatible with',('RFC10','RFC21'))
        enzyme_set=get_standard(standards)
        st.sidebar.markdown('')
        st.sidebar.subheader('Promoter')
        promoter_strength=st.sidebar.slider('Desired promoter strength (relative to BBa_J23100)',0.0,1.0,key='promoter_strength')
        st.sidebar.markdown('')
        st.sidebar.subheader('RBS')
        RBS_strength=st.sidebar.slider('Desired RBS strength (relative to BBa_B0034)',0.0,1.0,key='RBS_strength')
        st.sidebar.markdown('')
        st.sidebar.subheader('Target Sequence')
        target_sequence=st.sidebar.text_area('Enter your target sequence')
        if not Check_DNA(target_sequence):
                st.sidebar.markdown('Please enter a valid target sequence')
        st.sidebar.markdown('Enter the size of the pool to generate the toeholdswitch')
        pool=st.sidebar.number_input('Size of the pool to generate',3,10000,step=10)
        st.sidebar.markdown('Enter the minimum acceptable On/Off ratio')
        treshold=st.sidebar.number_input('Minimum acceptable On/Off ratio',10.0,70.0,step=5.0)
        st.sidebar.markdown('')
        st.sidebar.subheader('Output fron the toholdswitch')
        output_type=st.sidebar.selectbox('Type of output',('Reported BioBrick reporter protein','Your own output protein'))
        if output_type=='Your own output protein':
                output_sequence=st.sidebar.text_area('Enter your output sequence')
                if not Check_DNA(output_sequence):
                        st.sidebar.markdown('Please enter a valid output sequence')
                elif not test_standard(output_sequence,enzyme_set):
                        st.sidebar.markdown('Sequence not compatible with the selected standards')
                else:
                        output_df=pd.DataFrame({'Id':'NA','Description':'User-defined output sequence','Sequence':output_sequence},index=[0])
        else:
                reporter_df=load_reporters(standards)
                reporters=tuple(reporter_df['Description'])
                st.sidebar.markdown('Select the desired reporter protein')
                output_name=st.sidebar.selectbox('Reporter protein',reporters)
                reporter_index=reporter_df.index[reporter_df['Description']==output_name].tolist()
                reporter_index=int(reporter_index[0])
                output_df=reporter_df[reporter_index:reporter_index+1]
                output_df.columns=['Id','Description','Sequence','Standard']
        if st.sidebar.button('Assemble BioBrick',key='Assemble BioBrick'):
                promoter=Promoter_Selection(promoter_strength)
                rbs=RBS_Selection(RBS_strength)
                toehold=ToeholdSequence_gen(target_sequence,pool,rbs,treshold)
                final_seq,result_df=Assembler(promoter,rbs,toehold,output_df)
                st.subheader('Results')
                st.markdown('Plain Sequence')
                st.markdown('')
                st.write("5'-"+final_seq[:int(len(final_seq)/2)])
                st.write(final_seq[int(len(final_seq)/2):]+"-3'")
                st.markdown('')
                st.markdown('Annotation Table')
                st.write(result_df.set_index('Type'))
                st.markdown(get_table_download_link(result_df), unsafe_allow_html=True)
                st.markdown('Map')
                colors=['#ADD8E6','#00FF00','#00FFFF','#FFFFFF','#008000','#0000FF','#FFA500','#800080','#FF0000','#ADD8E6']
                features=[]
                for ind in list(result_df.index.values):
                        features.append(GraphicFeature(start=result_df.loc[ind,'Start'], end=result_df.loc[ind,'End'], strand=+1,color=colors[ind],label=result_df.loc[ind,'Type']))
                record=GraphicRecord(sequence_length=len(final_seq), features=features)
                record.plot(figure_width=5)
                st.pyplot()
                
        

        

if  __name__=='__main__':
        main()
