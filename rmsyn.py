import os

def main(f, L, D, I=None,Q=None,U=None,V=None, outp=''):
    d = ['P_rm', 'I_rm', 'V_rm']
    for i in range(3):
        try:
            os.system('rm -r %s/%s'%(outp,d[i]))
        except: None
    if Q != None and U != None:
        os.system('mkdir %s/%s'%(outp,d[0]))
        comm = 'rmsynthesis -o %s/%s --low -%s --high %s --dphi %s -m 11 %s %s %s' % \
           (outp,d[0], L,L,D, Q, U, f)
        print comm
        os.system(comm)
    elif I != None:
        os.system('mkdir %s/%s'%(outp,d[1]))
        comm = 'rmsynthesis -o %s/%s --low -%s --high %s --dphi %s -m 11 -q 1. -u 0. %s %s %s' % \
           (outp,d[1], L,L,D, I, I, f)
        os.system(comm)
    elif V != None:
        os.system('mkdir %s/%s'%(outp,d[2]))
        comm = 'rmsynthesis -o %s/%s --low -%s --high %s --dphi %s -m 11 -q 1. -u 0. %s %s %s' % \
           (outp,d[2], L,L,D, V, V, f)
        os.system(comm)

if __name__=='__main__':
    main()
