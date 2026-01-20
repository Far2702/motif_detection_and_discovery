from flask import Flask, request,redirect,render_template,url_for
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import logomaker

app=Flask(__name__)

# Homepage
@app.route('/',methods=['GET'])
def home():
    return render_template('home.html')

# Scanpage
@app.route('/result',methods=['POST'])
def result():
    
    sequence=request.form['sequence']
    alpha=request.form['alpha']
    alpha=float(alpha)
    
    # PPM I got from my previous calculations
    PPM=np.array([[0.44458545, 0.46573604, 0.38113367, 0.11886633, 0.07656514,
        0.08079526, 0.03426396, 0.77030457, 0.2034687 , 0.25      ,
        0.22038917, 0.40651438, 0.25423012, 0.44035533, 0.07233503,
        0.10617597, 0.71108291, 0.15693739, 0.69839255, 0.2965313 ,
        0.3642132 , 0.35152284],
       [0.40651438, 0.37267343, 0.33460237, 0.67301184, 0.17385787,
        0.70685279, 0.11463621, 0.04695431, 0.38536379, 0.24576988,
        0.30922166, 0.18654822, 0.19077834, 0.24153976, 0.6857022 ,
        0.11463621, 0.11463621, 0.07233503, 0.0892555 , 0.40651438,
        0.49111675, 0.48265651],
       [0.11463621, 0.12732657, 0.15270728, 0.07656514, 0.69416244,
        0.07656514, 0.81683587, 0.09348562, 0.17385787, 0.21615905,
        0.16962775, 0.25      , 0.35152284, 0.11463621, 0.11463621,
        0.04695431, 0.12732657, 0.04695431, 0.12732657, 0.14847716,
        0.08079526, 0.05541455],
       [0.03426396, 0.03426396, 0.13155668, 0.13155668, 0.05541455,
        0.1357868 , 0.03426396, 0.0892555 , 0.23730964, 0.28807107,
        0.30076142, 0.15693739, 0.2034687 , 0.2034687 , 0.12732657,
        0.7322335 , 0.04695431, 0.72377327, 0.08502538, 0.14847716,
        0.06387479, 0.11040609]])
    
    
    # Empirical Distirbtuion of Log-Likelihood Ratio Scores
    sampeled_seq=random_background_seq(10000)  # Number of smaples to be used is fixed by me

    log_like_ratio_score=[]

    for i in sampeled_seq:
        log_likelihood=0
        for j in range(len(i)):
            if i[j]=="A":
                log_likelihood=log_likelihood+np.log(PPM[0,j]/0.25)
            elif i[j]=="T":
                log_likelihood=log_likelihood+np.log(PPM[1,j]/0.25)
            elif i[j]=="G":
                log_likelihood=log_likelihood+np.log(PPM[2,j]/0.25)
            elif i[j]=="C":
                log_likelihood=log_likelihood+np.log(PPM[3,j]/0.25)

        log_like_ratio_score.append(log_likelihood)
    
    threshold=np.percentile(log_like_ratio_score,100*(1-alpha))
    
    
    

    
    scores={}
    i=0
    j=21
    
    while j<=len(sequence)-1:
        log_likelihood=0
        index=0
        
        for k in range(i,j+1):
            if sequence[k]=="A":
                log_likelihood=log_likelihood+np.log(PPM[0,index]/0.25)
            elif sequence[k]=="T":
                log_likelihood=log_likelihood+np.log(PPM[1,index]/0.25)
            elif sequence[k]=="G":
                log_likelihood=log_likelihood+np.log(PPM[2,index]/0.25)
            elif sequence[k]=="C":
                log_likelihood=log_likelihood+np.log(PPM[3,index]/0.25)
            index+=1
        
        scores[i]=log_likelihood
        i+=1
        j+=1
    
    
    significant_sites=[]
    for i in scores:
        if scores[i]>threshold:
            significant_sites.append(sequence[i:i+22])
    
    print("Significant Sites:",significant_sites)
    


    ppm_df=pd.DataFrame(
        PPM.T,               
        columns=['A','T','G','C']
    )
    
    # Logomaker plot
    logo=logomaker.Logo(
        ppm_df,
        shade_below=.5,
        fade_below=.5
    )

    logo.ax.set_xticks(range(len(ppm_df)))
    logo.ax.set_xticklabels(range(len(ppm_df)))

    logo.ax.set_xlabel("Position")
    logo.ax.set_ylabel("Probability")
    plt.title("Sequence Logo (CRP Motif)")
    plt.savefig("static/sequence_logo.png")
    plt.close()
    
    
    # Epmirical Distribution Plot
    plt.hist(log_like_ratio_score,edgecolor='black',density=True,bins=30)
    plt.title('Empirical Disitrbtuion of Log-Likelihood Ratio Score')
    plt.xlabel('S(X)')
    plt.ylabel('Density')
    plt.axvline(x=threshold,color='black',linestyle='--',label=f'threshold ({threshold})')
    plt.legend()
    plt.savefig("static/empirical_dist.png")
    plt.close()
    
    
    return render_template('results.html',significant_sites=significant_sites,alpha=alpha)
    

        
def random_background_seq(n):

    created_seq=[]

    for i in range(n):
        x=np.random.randint(0,4,size=(22,))
        seq=["a"]*22
        for j in range(len(x)):
            if x[j]==0:
                seq[j]="A"
            elif x[j]==1:
                seq[j]="T"
            elif x[j]==2:
                seq[j]="G"
            elif x[j]==3:
                seq[j]="C"

        created_seq.append(seq)

    return created_seq



if __name__=='__main__':
    app.run(debug=True)
    