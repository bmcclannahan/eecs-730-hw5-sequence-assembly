from align_sequence import local_alignment

class Node:
    def __init__(self,sequence):
        self.sequence = sequence
        self.next = None
        self.next_score = 0
        self.previous = None
        self.previous_score = 0

    def set_next(self,node,score):
        self.next = node
        self.next_score = score
    
    def set_previous(self,node,score):
        self.previous = node
        self.previous_score = score

def get_start_node(nodes):
    m = 10000
    start = None
    for node in nodes:
        #print(node.sequence,node.prev_score)
        if node.previous_score < m:
            start = node
            m = node.previous_score
    return start

def get_end_node(nodes):
    m = 10000        
    end = None
    for node in nodes:
        #print(node.sequence,node.next_score)
        if node.next_score < m:
            end = node
            m = node.next_score
    return end

def print_graph(start,end):
    curr = start
    while curr.sequence != end.sequence:
        print(curr.sequence)
        curr = curr.next
    print(end.sequence)

def build_sequence(start,end):
    while start.next.sequence != end.sequence:
        start = merge(start,start.next)
    complete = merge(start,end)
    return complete

def get_max_index(matrix):
    curr = matrix[0][0]
    index = (0,0)
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > curr:
                index = (i,j)
                curr = matrix[i][j]
    return index

def build_graph(sequences,necessary_score = 0):
    nodes = []
    for sequence in sequences:
        node = Node(sequence)
        nodes.append(node)
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i == j:
                continue
            mat = local_alignment(nodes[i].sequence,nodes[j].sequence,mismatch_cost=-1)
            f,l = get_max_index(mat)
            m = mat[f][l]
            if f > l and m >= necessary_score and nodes[i].next_score < m:
                nodes[i].set_next(nodes[j],m)
                nodes[j].set_previous(nodes[i],m)
            elif l < f and m >= necessary_score and nodes[i].previous_score < m:
                nodes[j].set_next(nodes[i],m)
                nodes[i].set_previous(nodes[j],m)
    return nodes

def merge(node1,node2):
    new_seq = ''
    if len(node1.sequence) == len(node2.sequence):
        i = 0
        while node1.sequence[i:] != node2.sequence[:0-i]:
            i+=1
        new_seq = node1.sequence[0:i] + node2.sequence
    elif len(node1.sequence) > len(node2.sequence):
        a = len(node1.sequence) - len(node2.sequence)
        i = 0
        while node1.sequence[a+i:] != node2.sequence[:0-i]:
            i+=1
        new_seq = node1.sequence[0:i+a] + node2.sequence
    else:
        a = len(node2.sequence) - len(node1.sequence)
        i = 0
        while node1.sequence[a+i:] != node2.sequence[:0-i-a]:
            i+=1
        new_seq = node1.sequence[0:i] + node2.sequence
    node = Node(new_seq)
    node.set_previous(node1.previous,node1.previous_score)
    node.set_next(node2.next,node2.next_score)
    return node