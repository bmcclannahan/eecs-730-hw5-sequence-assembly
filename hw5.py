from align_sequence import local_alignment, print_alignment_matrix
from sequence_assembler import build_graph, merge, get_start_node, get_end_node, print_graph, build_sequence
from read_fasta import read_file
sequence = 'atcgatcattact'
sequences = [sequence[0:6],sequence[2:8],sequence[4:10],sequence[1:7],sequence[7:],sequence[6:12]]
print(sequences)

reads = read_file('HW5_reads.fasta')
sequences = reads

nodes = build_graph(sequences,necessary_score=50)

start = get_start_node(nodes)
end = get_end_node(nodes)

f = open('hw5.out','a+')
f.write(start.sequence+'\n')
f.write(end.sequence+'\n')
curr = start
while curr.sequence != end.sequence:
    f.write(curr.sequence+'\n')
    curr = curr.next
print(end.sequence+'\n')

#print(start.sequence,end.sequence)
#print_graph(start,end)

complete = build_sequence(start,end)
#print(complete.sequence)

f.write('\n' + complete.sequence)