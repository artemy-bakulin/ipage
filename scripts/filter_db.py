import MI
import numpy as np

def non_redundancy_sort_pre(go_profiles):
    go_sums={go:go_profiles[go].sum() for go in go_profiles} 
    go_profiles={go:go_profiles[go] for go in go_profiles if go_sums[go]>6}
    syn_ungrouped=[]
    f1=lambda go2: go1!=go2 and  x1<go_sums[go2]<x2 
    f2=lambda go2: np.count_nonzero((go_profiles[go1]+go_profiles[go2]) == 2)>min(go_sums[go1],go_sums[go2])*0.6
    f3=lambda go2: MI.mut_info_normalized(go_profiles[go1],go_profiles[go2],2,2)>0.1
    go_li=set()
    groups=[]
    for go1 in go_profiles:
        if go1 not in go_li:
            x1=0.5*go_sums[go1]
            x2=1.5*go_sums[go1]
            group={go2 for go2 in {go2 for go2 in {go2 for go2 in go_profiles if f1(go2)} if f2(go2)} if f3(go2)}
            groups.append(group)
            go_li=go_li|group
    groups=set([i for el in groups for i in el ])
    for el in groups:
        del go_profiles[el] 
    return go_profiles
    
def non_redundancy_sort_post(cmi_output, db_profiles):
    mi_dic={(name1,name2):MI.mut_info(db_profiles[name1],db_profiles[name2],2,2) for name1 in cmi_output for name2 in cmi_output}
    sorted_mi=sorted(mi_dic,key=lambda x: mi_dic[x],reverse=1)
    max_index=sorted_mi.index(list(filter(lambda x: x[0]==x[1], sorted_mi))[-1])+1
    connections=list(filter(lambda x: x[0]!=x[1] and mi_dic[x]>0.025,sorted_mi))
    names={el[0] for el in mi_dic}
    i=0
    for el in connections:
        if i%2==0:
            if el[0] not in names or el[1] not in names or abs(sum(db_profiles[el[0]])-sum(db_profiles[el[1]]))<50:
                pass
            elif sum(db_profiles[el[0]])<sum(db_profiles[el[1]]):
                names.remove(el[1])
            else:
                names.remove(el[0])
        i+=1
    names=names|({el[0] for el in mi_dic}-{el[0] for el in connections})
    return {el:cmi_output[el] for el in names}