import pywt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import sys
import os
def get_range(df, c, comparison):
    datac = df[df['r'] < c if comparison == '<' else df['r'] > c].index.tolist()
    if not datac:
        return []

    result = [[datac[0]]]
    for i in range(1, len(datac)):
        if datac[i - 1] + 1 == datac[i]:
            result[-1].append(datac[i])
        else:
            result.append([datac[i]])

    fin_range = [[df.loc[i[0], 'loc'], df.loc[i[-1], 'loc']] for i in result]
    return fin_range

def get_crange(df, c=-.4):
    return get_range(df, c, '<')

def get_hrange(df, c=0.6):
    return get_range(df, c, '>')


def wavelet_transform(df, rlist,a,hc=0.6,lc=-0.4):
    rec=[float(x) for x in rlist]
    if a == "30":
        wavefunc = 'db4'
        level = 5
        threshold = 0.8
        
    else :
        wavefunc= 'db2'
        level = 8
        threshold = 0.8
    coeffs = pywt.wavedec(rec, wavefunc,mode='periodic', level=level)
    for i in range(1, len(coeffs)):
        coeffs[i] = pywt.threshold(coeffs[i], threshold * max(coeffs[i]))     
    y_denoised = pywt.waverec(coeffs, wavefunc,mode='periodic',)
    if df.shape[0] %2 !=0:
        y_denoised=y_denoised[:-1]
    y_re=pd.DataFrame(y_denoised,columns=['r'],index=range(df.shape[0]))
    y_re['loc']=df.index
    if a == '30':
        coldspot=y_re[y_re['r']>hc]['loc'].tolist()
    else:
        coldspot=y_re[y_re['r']<lc]['loc'].tolist()
    return rec,y_denoised,y_re,coldspot

def plot_and_save(get_crange, get_hrange, wavelet_transform, result_path, widthf=6/2.54,hc=0.6,lc=-0.4):
    # 导入数据集
    params = {
            "font.size": 6,     # 全局字号
            "axes.spines.right":False,  # 坐标系-右侧线
            "axes.spines.top":False,   # 坐标系-上侧线
            "axes.titlesize":10,   # 坐标系-标题-字号
            "axes.labelsize": 8,  # 坐标系-标签-字号
            "legend.loc":'upper right',  # 图例-位置
            "legend.fontsize": 8,  # 图例-字号
            "xtick.labelsize": 6,  # 刻度-标签-字号
            "ytick.labelsize": 6,  # 刻度-标签-字号
            "xtick.direction":'in',   # 刻度-方向
            "ytick.direction":'in',  # 刻度-方向
            'pdf.fonttype':42,



    }
    mpl.rcParams.update(params)
    for a in ['30','all']:
        df=pd.read_table(result_path,index_col=0,header=0)
        with open(f'{result_path}.fa.pmr.{a}.2.txt','r') as f:
            rlist=f.read().rstrip().splitlines()
        rec, y_denoised ,y_re,coldspot = wavelet_transform(df, rlist,a,hc,lc)
        plt.figure().set_size_inches(widthf,widthf*0.618)
        plt.ylim(-1, 1)
        plt.tick_params(which='major',width=0.5,length=3)

        with open(f'{result_path}.fa.pmr.{a}.coldsnp.txt','w')as f:
            f.write('\n'.join([str(x) for x in coldspot]))
        plt.scatter(df.index.tolist(), rec, label='original signal',s=0.3,c="none",alpha=0.1,marker='o',edgecolors="grey",linewidths=0.3)
        plt.plot(df.index.tolist(), y_denoised, label='Denoised signal',linewidth=0.5)
        plt.xlabel('Genome Position (Mb)')
        plt.ylabel('R value')
        if a=='30':
            hotspot=get_hrange(y_re,hc)
            n=0
            for i in hotspot:
                if n==0:
                    n+=1
                    plt.fill_between([i[0],i[1]],0.95,1,facecolor = '#FF0000' ,label='hot spots')
                    continue
                plt.fill_between([i[0],i[1]],0.95,1,facecolor = '#FF0000', )
            hotspot=y_re[y_re['r']>hc]['loc'].tolist()
            with open(f'{result_path}.fa.pmr.{a}.hotsnp.txt','w')as f:
                f.write('\n'.join([str(x) for x in hotspot]))

        if a=="all":
            n=0
            coldspot=get_crange(y_re,lc)
            for i in coldspot:
                if n==0:
                    n+=1
                    plt.fill_between([i[0],i[1]],-1,-0.95,facecolor = '#0000FF' ,label='cold spots')
                    continue
                plt.fill_between([i[0],i[1]],-1,-0.95,facecolor = '#0000FF', )
        if a=='30':
            plt.axhline(y=hc,c='r',linewidth=0.5)
        else:
            plt.axhline(y=lc,c='r',linewidth=0.5)
        plt.savefig(f'{result_path}.fa.pmr.{a}.2.wt.l.pdf',dpi=300,bbox_inches='tight')
if __name__ == "__main__":
    result_path = sys.argv[1]
    h_criterion=float(sys.argv[2])
    l_criterion=float(sys.argv[3])
    if not os.path.exists(result_path):
        raise FileNotFoundError(f"File '{result_path}' not found.")
    widthf=6/2.54
    plot_and_save(get_crange, get_hrange,  wavelet_transform,result_path=result_path,hc=h_criterion,lc=l_criterion)
