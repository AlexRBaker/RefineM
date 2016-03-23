###############################################################################
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import matplotlib
import mpld3
import logging

import numpy as np
import os
import sys
import biolib.seq_io as seq_io

#Changed the direciton of this library. Now, it is going to use biopython
#to read in the fast afile, make windows and then write a new fasta file as well as a links file
#these are piped into refinem

class WindowGen(object):
    def __init__(self,cpus):
        self.cpus=cpus
        
    def write_windows(self,scaffold_file,output_dir,window_size,window_gap):
        '''
        --------------------------------------------------------------------
        Take a scaffold file in fasta format and prints a similarly name
        fasta file of the windows made from the scaffolds in the scaffold file.
        --------------------------------------------------------------------
        Input: scaffold_file
                    The name of the fasta file to turn into windows
        Output:
                Writes a links file and a file of windows
        '''
        seq_win_id={} #pairs a scaffold with the windows made from it
        window_dict={} #dictionary of windows_dict[win_id]=seq_win
        
        for seq_id, sequence in seq_io.read_seq(scaffold_file):
            win_id,seq_win=self.make_windows([seq_id,sequence], window_size, window_gap)
            seq_win_id[seq_id]=win_id
            for i in range(0,len(win_id)):
                window_dict[win_id[i]]=seq_win[i]
        
        
        filename=os.path.split(scaffold_file)[1]
        start="".join(filename.split('.')[:-1])
        end=filename.split('.')[-1]
        window_file=os.path.join(output_dir,start+"windows."+end)    
        
        names = ['id','sequence']
        formats = []
        
        
        print len(window_dict)
        seq_io.write_fasta(window_dict,window_file)
        #self.write_fasta(window_dict,window_file)
            
        links_file=os.path.join(output_dir,"links_file.tsv")
        self.write_links(seq_win_id,links_file)
        return [window_file,links_file]
        
    #~ def write_fasta(self,dictionary,file_name):
        #~ 
        #~ #memory intensive
        #~ fasta=np.array([">{0}{1}{2}{3}".format(fasta_id,os.linesep,sequence,os.linesep) for fasta_id,sequence in dictionary.iteritems()])
        #~ np.savetxt(file_name,fasta,newline='')
            #~ #fasta="".join(fasta)
            #~ #fasta_file.write(fasta)
#~ 
            #~ for fasta_id, sequence in dictionary.iteritems():
                #~ fasta_file.write(">{0}{1}{2}{3}".format(fasta_id,os.linesep,sequence,os.linesep))
        #~ #not memory intensive
        #~ with open(file_name) as fasta_file:
            #~ for fasta_id, sequence in dictionary.items():
                #~ fasta_file.write(">{0}{1}{2}{3}".format(fasta_id,os.linesep,sequence,os.linesep))
         
    def make_windows(self,seq_info,window_size,window_gap):
        '''
    -------------------------------------------------------------
    Makes a series of sliding windows from a sequence based on
    the chosen windows size and gap distance
    --------------------------------------------------------------
    Input:   seq_info
                seq_id
                    The unique ID for that sequence
                sequence
                    The sequence in nucleotide space
             window_size
                The size of the window to slide across the sequence
                it has two modes based on the type of the input
                    Int: Make windows that many characters long
                    float: Make windows equal to that proportion of
                        the total length
             gap_size
               The size of the gap between windows
                    Int: Gaps will be that many chracters long
                    float: Gaps will be that proportion of the total 
                        sequence length
    Output:
        Windows
            win_id
                A unique modified seq_id for each window
            windows
                A list of all of the windows produced by the specified 
                windows size and gap. If the last window is larger than
                the remaining sequence then the all remaining sequence
                is used.
        '''
        seq_id,sequence=seq_info
        seq_length=len(sequence)
        win_seq_id={}
        
        if isinstance(window_size,str):
            window_size=self.type_convert(window_size)
        elif not isinstance(window_size,int) and not isinstance(window_size,float):
            raise TypeError
            
        if isinstance(window_gap,str):
            window_gap=self.type_convert(window_gap)
        elif not isinstance(window_gap,int) and not isinstance(window_gap,float):
            raise TypeError
            
        try:
            if isinstance(window_size,int) and isinstance(window_gap,int):
                #No changed to be made - non-parametrised step, gap assumed
                pass
            elif isinstance(window_size,int) and isinstance(window_gap,float):
                #If not an int then fractional of total length is assumed
                window_gap=int(window_gap*seq_length)
            elif isinstance(window_size,float) and isinstance(window_gap,int):
                #Assumed the window_size is a fraction of total length
                window_size=int(window_size*seq_length)
            elif isinstance(window_size,float) and isinstance(window_gap,float):
                #Both are assumed to be fraction of total sequence length
                window_size=int(window_size*seq_length)
                window_gap=int(window_gap*seq_length)
            else:
                raise TypeError
            Windows=[[],[]]
            for i in range(0,seq_length,window_size+window_gap):
                if len(sequence[i:(i+window_size)])>=100:
                    Windows[0].append("{0}:{1}to{2}".format(seq_id,i,min(i+window_size,seq_length)))
                    Windows[1].append(sequence[i:(i+window_size)])
            return Windows #returns [list of ids, list of windows]
        except TypeError:
            print "The window size and gap must be int(#bp) or float(proportion of total length)"
            raise
            
    def type_convert(self,num_str):
        '''Covert a number string into appropiate int or float '''
        try:
            return int(num_str)
        except ValueError:
            return float(num_str)
        except TypeError:
            print "Inputted object type:",type(num_str)
            raise
        
    def write_links(self, Windows,links_file):
        '''
        ----------------------------------------------------------
        Creates a links file for use by refinem. The links will be
        between the windows made from each contig. The file is just 
        tab separated id names telling refinem what to link.
        -----------------------------------------------------------
        
        Input:
            Windows: dict[scaffold_name]=[list of created windows ids]
                scaffold_name
                    The name of the scaffold used to make the windows
                windows_ids
                    The id of windows created from this scaffold
                    
        Output:
            Boolean
                True if successfully wrote the file
                otherwise raises an error
        '''
        complete=False
        try:
            index=range(0,len(Windows))
            with self.tryopen(links_file) as linksfile:
                for scaffold,windows in Windows.iteritems():
                    for i in range(0, len(windows)-1):
                        linksfile.write("{0}\t{1}\n".format(windows[i],windows[i+1]))
                complete=True
            return complete
        except:
            raise

    def tryopen(self, filename):
        '''
        -----------------------------------------------------------
        Opens a file. If it does not exist then it creates a file
        with the specified name.
        -----------------------------------------------------------
        
        Input: 
            filename
                name of the file
        Output:
            opened file
        '''
        try:
            if not os.path.isfile(filename):
                with open(filename,'w') as window_links:
                    pass
            else:
                print "This file already exists. New content will be appended"
            return open(filename,'a+')
        except IOError:
            print "There was an error in reading or writing the file"
            raise
        except:
            print "There was an unidentified error"

    #~ def run(self, seq_info,window_size,window_gap,filename):
        #~ '''
    #~ --------------------------------------------------------------------------------------- 
        #~ runs the sliding windows, writes the links and returns a set of windows with
        #~ individual ids.
    #~ ---------------------------------------------------------------------------------------
    #~ Output:
        #~ window_info
            #~ win_id
                #~ Unique id associated with the window
            #~ window
                #~ The window of the overall sequence taken
                #~ 
        #~ Also, a links_file is written with user specified name.
        #~ '''
        #~ new_windows=windows(seq_info, window_size,window_gap)
        #~ 
        #~ write_links(new_windows,filename)
        #~ 
        #~ return new_windows #lists of ids and windows [[ids],[windows]]
    
