import sys


def cigarToList(cigar,seq):
    ret, i = [], 0
    op_map = {'M':0, # match or mismatch
              '=':0, # match
              'X':0, # mismatch
              'I':1, # insertion in read w/r/t reference
              'D':2, # deletion in read w/r/t reference
              'N':3, # long gap due e.g. to splice junction
              'S':4, # soft clipping due e.g. to local alignment
              'H':5, # hard clipping
              'P':6} # padding
  
    while i < len(cigar):
        run = 0
        while i < len(cigar) and cigar[i].isdigit():
            # parse one more digit of run length
            run *= 10
            run += int(cigar[i])
            i += 1
        assert i < len(cigar)
        # parse cigar operation
        op = cigar[i]
        i += 1
        #print("op=",op)
        #print("op_map=",op_map)
        assert op in op_map
        # append to result
        ret.append([op_map[op], run])
    new_seq=gapped_seq(ret,seq)    
    return ret,new_seq




def gapped_seq(new_cgar,seq):
        edit_seq=''
        remain=''
        new_seq=''
        del_seq=''
        remain=seq

        for a,b in new_cgar:
                if a==1:  #Insertion
                        extract = len(remain) - b
                        tmp = remain[-extract:]
                        remain = tmp
                        continue

                if a==2:  #Deletion
                        del_seq = int(b) * '-'
                        edit_seq = edit_seq + del_seq
                        continue

                else:
                        new_seq = remain[:b]
                        edit_seq = edit_seq + new_seq
                        chopping_len = len(remain)-b
                        remain = seq[-chopping_len:]
                        #print("new", edit_seq)
                        #print("remain",remain)

        return edit_seq




#cigar = sys.argv[1] 
#seq = sys.argv[2]
#new_cgar,new_seq=cigarToList(cigar,seq)#print(cigar)
#print(new_cgar)
#print(seq)
#print(new_seq)

