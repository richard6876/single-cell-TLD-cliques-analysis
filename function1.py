import copy

#read TLD results file of deDoc2.w
def read_TLDfile_de2w(chr_name,level):
    '''
    :param chr_name:  chromosome name
    :param level:  sparsity level 
    :return ls:  TLD list
    '''
    chr_namename="chr_name"+chr_name
    f=open("./deDoc2s/RaoIMR90CallTAD/"+level+"/raw/binSize40000/"+chr_namename+"_10kb.RAWobserved."+level+".matrix.window.TAD","r")
    ls=[]
    for line in f:
        x=line.strip('\n').split(' ')
        #print(len(x))
        if len(x)!=2:
            x1=40000*(eval(x[0])-1)
            x2=40000*eval(x[-2])
            ls.append([chr_namename,str(x1),str(x2)])
    return ls

#read contact file
def read_contact_file(chr_name,level):
    '''
    :param chr_name:  chromosome name
    :param level:  sparsity level 
    :return ls:  contact list
    '''
    chr_namename="chr_name"+chr_name
    f=open("/home/wanghy/deDoc2_clique/RaoIMR90/RAWobserved-"+level+"/"+chr_namename+"_10kb.RAWobserved."+level,"r")
    ls=[]
    for line in f:
        ls.append(line.strip('\n').split('\t'))
    f.close()
    return ls

#get contact file per chromosome
def get_chr_contact(chr_name,contact_list):
    '''
    :param chr_name:  chromosome name
    :param contact_list:  contact list of all chromosomes 
    :return ls_c:  contact list of the specific chromosome
    '''
    ls_c=[]
    for row in contact_list:    
        if row[0]==chr_name and row[2]==chr_name:
            ls_c.append(row)
    return ls_c
    
# calculate the TLD-TLD contacts
def cal_contact1(TLDlist,contactlist):
    '''
    :param TLDlist:  TLD coordinate list
    :param contactlist:  contact list of the specific chromosome 
    :return TLD_contacts_list:  contact list of the TLD-TLD on the specific chromosome
    '''  
    TLD_contacts_list=[]
    
    for contact in contactlist:
        tad1=[]
        tad2=[]
        for tad in TLDlist:
            if eval(tad[1])<=eval(contact[0]) and eval(tad[2])>=eval(contact[0]):
                tad1=tad
                for t in TLDlist:
                    if eval(t[1])<=eval(contact[1]) and eval(t[2])>=eval(contact[1]):
                        tad2=t
                if tad1 != [] and tad2 != [] and tad1 != tad2:
                    TLD_contacts_list.append([tad1,tad2,contact[2]])
                    
    TLD_contacts_list_1=copy.deepcopy(TLD_contacts_list)
    for i in TLD_contacts_list:
        i[2]=str(eval(i[2])-1)
        for j in TLD_contacts_list_1:
            if (i[0]==j[0] and i[1]==j[1]) or (i[0]==j[1] and i[1]==j[0]):
                if j[2]=='10':
                    TLD_contacts_list.remove(i)
                    break
                else:
                    i[2]=str(eval(i[2])+eval(j[2]))
                    j[2]='10'  
    return TLD_contacts_list

#write the information to bedfile
def write_bedfile(software,TLD_contacts_list,level):
    '''
    :param software:  TLD coordinate list
    :param TLD_contacts_list:  contact list of the TLD-TLD on the specific chromosome
    :param level:  sprsity level of the data
    :return :  bedfile
    '''  
    f=open("/home/wanghy/deDoc2_clique/results/"+software+"."+level+".TAD_interaction.bed","w")
    for row in TLD_contacts_list:
        x=[]
        x.extend(row[0])
        x.extend(row[1])
        x.append(str(int(eval(row[2]))))
        f.write("\t".join(x)+"\n")
    f.close()
