export = open('config.div','w')
conf_file = open('config.div','r')
data = conf_file.readlines()
conf_file.close()

directory = '../../TrainingHalos/data/particle_ID/MDmini'
x,y,z =  '2','3','4'
skip_lines = '16'
processes = '4'
steps = '50000'
plot = 'False'
test = 'False'

line =  directory+' '+x+' '+y+' '+z+' '+skip_lines+' '+processes+' '+steps+' '+plot+' '+test+'\n'

export.write(line)
export.close()