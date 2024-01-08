import time
import copy
from function1 import read_celllist,read_TLDfile_de2w,read_contact_file,get_chr_contact,cal_contact1,write_bedfile

# get the TLD-TLD contact list per cell
if __name__ == '__main__':
    start = time.time()
    cell_list=read_celllist()
    chr_name_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
    software="de2w"
    for cell in cell_list:
        contact_list=read_contact_file(cell)
        f_list=[]
        for chr_name in chr_name_list:
            try:
                contact_chr_list=get_chr_contact(chr_name,contact_list)
                TADlist=read_TLDfile_de2w(cell,chr_name)
            except:
                continue
            clist=cal_contact1(TADlist,contact_chr_list)
            f_list=f_list+clist
        write_bedfile(software,f_list,cell)
    end = time.time()
    print(end-start)