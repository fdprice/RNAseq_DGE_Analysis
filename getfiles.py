import os
# ['M12_N706', '150723_NS500422_0164_AHF5C2BGXX'], ['M9_N705', '150722_NS500422_0163_AHFGFHBGXX'],
for prefix, project in [ ['M1_N701', '150806_NS500422_0172_AHCTCYBGXX'] ]:
    for i in [1, 2]:
        output_filename = '%s_run1_R%d.fastq.gz' % (prefix, i)
        url = 'http://software.rc.fas.harvard.edu/ngsdata/%s/All_Reads/Fastq/Undetermined_S0.R%d.fastq.gz' % (project, i)
        print url, output_filename
        os.system('curl -o %s %s' % (output_filename, url))


