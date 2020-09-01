import sys, getopt, commands

def main(argv):
   indir = ''
   outdir = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["indir=","outdir="])
   except getopt.GetoptError:
      print 'test.py -i <inputdir> [-o <outputdir>]'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--indir"):
         indir = arg
      elif opt in ("-o", "--outdir"):
         outdir = arg
   print 'Input Directory is "', indir
   print 'Output Directory is "', outdir
   samples=commands.getoutput('ls '+indir+'/*')
   print samples
   for i in range(0,len(samples)):
       #os.system('hadd ' + outdir + '/' +  +'.root ' + indirec +'/*' +'root')
       print i

if __name__ == "__main__":
   main(sys.argv[1:])
