import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
def sample_statistics(file_name):
    data_WF = pd.read_csv("input_files/stat/" + file_name + '.csv')
    data_WF = data_WF[1:] #starting from first now; remove zero th row that is total
    print(data_WF)
    plt.figure()
    plt.plot(data_WF['LOS_bin'], data_WF['Uni_Enct_N'], color='tab:blue', label='Encounter', marker='x', markersize=2,
             linewidth=1)  # graph edges number
    plt.yscale('log') # use log scale
    #plt.title('Number of edges of each ego network of C.diff')
    plt.legend()
    plt.xticks(np.arange(1, 27, 2)) #x axis labels
    #plt.subplot(212)  # for nodes
    plt.plot(data_WF['LOS_bin'], data_WF['Uni_Pat_N'],color='tab:orange',label='Patient',marker='x',markersize=2, linewidth=1)  # graph nodes number;
    plt.legend()
    plt.plot(data_WF['LOS_bin'], data_WF['Dia_N'],color='tab:green',label='Diagnose',marker='x',markersize=2, linewidth=1)  # graph nodes number;
    plt.legend()
    plt.plot(data_WF['LOS_bin'], data_WF['Avg_Pat_Age'],color='tab:red',label='Average Age',marker='x',markersize=2, linewidth=1)  # graph nodes number;
    plt.legend()
    plt.xlabel('Time window')

    #plt.yscale('linear')
    plt.subplots_adjust(top=0.95, bottom=0.09, left=0.06, right=0.98, hspace=0.25,
                        wspace=0.45)
    plt.savefig(os.path.join("output_files/figures", file_name+ ".eps"), format="eps",dpi =1000)
    plt.show()

def SCI_chart(filename):
    data_wf = pd.read_csv("input_files/comorbidity_All_SCI/comorbidity_5/" + filename + '.csv')
    #draw pie chart
    less0_05 = data_wf[(data_wf['SCI'] < 0.05)].shape[0]
    less0_1 = data_wf[(data_wf['SCI'] >= 0.05) & (data_wf['SCI'] <= 0.1)].shape[0]
    greater0_1 = data_wf[(data_wf['SCI'] > 0.1)].shape[0]
    my_data = [less0_05,less0_1,greater0_1]
    my_labels = 'SCI: [0, 0.05)', 'SCI: [0.05, 0.1]', 'SCI: [0.1, 1]'
    plt.pie(my_data, labels=my_labels, autopct='%1.1f%%')
    #plt.title('My Title')
    plt.axis('equal')
    plt.subplots_adjust(top=0.99, bottom=0.01, left=0.01, right=0.99, hspace=0.25,
                        wspace=0.45)
    plt.savefig(os.path.join("output_files/figures", filename+ "_pie.pdf"), format="pdf",dpi =1000)
    #plt.show()
    #draw histogram
    df_hist = data_wf[(data_wf['SCI'] >= 0.05) & (data_wf['SCI'] <= 0.1)][['SCI']]
    print(df_hist)
    df_hist.plot.hist(grid = False,legend=False, bins = 50, rwidth = 0.98, color = '#33B5FF')
    plt.xlabel('SCI')
    plt.ylabel('Frequency')
    plt.subplots_adjust(top=0.97, bottom=0.10, left=0.1, right=0.98, hspace=0.25,
                        wspace=0.45)
    plt.savefig(os.path.join("output_files/figures", filename + "_histogram.pdf"), format="pdf", dpi=1000)
    plt.show()




if __name__ == "__main__":
    SCI_chart('5_comorbidity_network') #file name without .csv extension