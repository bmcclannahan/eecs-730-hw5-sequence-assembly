#from HW1
def global_alignment(seq1, seq2, process = None):
    matrix = [[0]*(len(seq2)+1) for i in range(len(seq1)+1)]
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if i == 0:
                matrix[i][j] = j*-3
            elif j == 0:
                matrix[i][j] = i*-3
            else:
                d = -2
                if seq1[i-1] == seq2[j-1]:
                    d = 1
                left = matrix[i-1][j] - 3
                up = matrix[i][j-1] - 3
                diag = matrix[i-1][j-1] + d
                matrix[i][j] = max([left,up,diag])
    if process is not None:
        return process.memory_info()[0]/1000000
    return matrix

def affine_gap_alignment(seq1, seq2):
    top_matrix = [[0]*(len(seq2)+1) for i in range(len(seq1)+1)]
    mid_matrix = [[0]*(len(seq2)+1) for i in range(len(seq1)+1)]
    bot_matrix = [[0]*(len(seq2)+1) for i in range(len(seq1)+1)]
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if i == 0 and j == 0:
                top_matrix[i][j] = -5
                mid_matrix[i][j] = 0
                bot_matrix[i][j] = -5
            elif i == 0:
                bot_matrix[i][j] = max(bot_matrix[i][j-1]-1,mid_matrix[i][j-1]-6)
                top_matrix[i][j] = (j)*-1-5
                mid_matrix[i][j] = max([top_matrix[i][j],top_matrix[i][j]])
            elif j == 0:
                bot_matrix[i][j] = (i)*-1-5
                top_matrix[i][j] = max(top_matrix[i-1][j]-1,mid_matrix[i-1][j]-6)
                mid_matrix[i][j] = max([bot_matrix[i][j],bot_matrix[i][j]])
            elif i != 0 and j != 0:
                top_matrix[i][j] = max(top_matrix[i-1][j]-1,mid_matrix[i-1][j]-6)
                bot_matrix[i][j] = max(bot_matrix[i][j-1]-1,mid_matrix[i][j-1]-6)
                match = -2
                if seq1[i-1]==seq2[j-1]:
                    match = 1
                mid_matrix[i][j] = max([mid_matrix[i-1][j-1]+match,top_matrix[i][j],bot_matrix[i][j]])
    return (top_matrix, mid_matrix, bot_matrix)

def local_alignment(seq1, seq2,mismatch_cost=-2):
    matrix = [[0]*(len(seq2)+1) for i in range(len(seq1)+1)]
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if i == 0:
                matrix[i][j] = 0
            elif j == 0:
                matrix[i][j] = 0
            else:
                d = mismatch_cost
                if seq1[i-1] == seq2[j-1]:
                    d = 1
                left = matrix[i-1][j] - 3
                up = matrix[i][j-1] - 3
                diag = matrix[i-1][j-1] + d
                matrix[i][j] = max([left,up,diag,0])
    return matrix

def indexed_alignment(query, reference, kmer = 12):
    kmer_dict = dict()
    for i in range(len(query)-kmer+1):
        k = query[i:i+kmer]
        if k in kmer_dict.keys():
            kmer_dict[k].append(i)
        else:
            kmer_dict[k] = [i]

    score_dict = dict()
    best_score = ([0],-5*len(query))
    for i in range(len(reference)-kmer+1):
        k = reference[i:i+kmer]
        if k in kmer_dict.keys():
            locs = kmer_dict[k]
            for loc in locs:
                start = i - loc
                if start < 0:
                    start = 0
                end = start + len(query)
                if end > len(reference):
                    end = len(reference)
                seq1 = query
                seq2 = reference[start:end]
                #print("Kmer match found with kmer:", k)
                #print("Sequence 1:", seq1)
                #print("Sequence 2:", seq2)
                top,mid,bot = affine_gap_alignment(seq1,seq2)
                score = mid[-1][-1]
                if best_score[1] < score:
                    best_score = ([start],score)
                elif best_score[1] == score and start not in best_score[0]:
                    best_score[0].append(start)
                score_dict[start] = mid[-1][-1]
                #print_affine_gap_alignment(top,mid,bot,seq1,seq2)

    # for key in score_dict.keys():
    #     print("Match found at ", key, " with alignment score of", score_dict[key])
    location,score = (best_score[0][0],best_score[1])
    print("The best match was found from", location+1, "to", location+len(query), "with a score of", score)
            
def print_global_alignment(alignment_matrix,seq1,seq2):
    first_sequence = ""
    second_sequence = ""
    gap_line = ""
    i = len(alignment_matrix)-1
    j = len(alignment_matrix[0])-1
    while i > 0 or j > 0:
        # print("i:",i,"j:",j)
        # print(first_sequence)
        # print(gap_line)
        # print(second_sequence)
        current = alignment_matrix[i][j]
        if i > 0 and current - alignment_matrix[i-1][j] == -3:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = "-" + second_sequence
            gap_line = " " + gap_line
            i -= 1
            # print(1)
        elif current - alignment_matrix[i][j-1] == -3:
            first_sequence = "-" + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = " " + gap_line
            j -= 1
            # print(2)
        elif i > 0 and j > 0 and current - alignment_matrix[i-1][j-1] == 1:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = "|" + gap_line
            j -= 1
            i -= 1
            # print(3)
        elif i > 0 and j > 0 and current - alignment_matrix[i-1][j-1] == -2:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = " " + gap_line
            j -= 1
            i -= 1
            # print(4)
        else:
            i -= 1
            j -= 1
            print(5)
    print("\nThe alignment is:")
    for i in range(len(first_sequence)//70+1):
        if len(first_sequence) < (i+1)*70:
            print(first_sequence[i*70:].rstrip())
            print(gap_line[i*70:].rstrip())
            print(second_sequence[i*70:].rstrip())
        else:
            print(first_sequence[i*70:(i+1)*70].rstrip())
            print(gap_line[i*70:(i+1)*70].rstrip())
            print(second_sequence[i*70:(i+1)*70].rstrip())

def print_affine_gap_alignment(top,mid,bot,seq1,seq2):
    current_matrix = 1 #0 is top, 1 is middle, 2 is bottom
    first_sequence = ""
    second_sequence = ""
    gap_line = ""
    i = len(mid)-1
    j = len(mid[0])-1
    # prev_cm = 1
    # prev_i = i
    # prev_j = j
    while i > 0 or j > 0:
        # if prev_cm != current_matrix or prev_i != i or prev_j != j:
        #     print(current_matrix,i,j)
        #     prev_cm = current_matrix
        #     prev_i = i
        #     prev_j = j
        if current_matrix == 0:
            #Extend Gap
            if top[i][j] == top[i-1][j] - 1:
                first_sequence = seq1[i-1] + first_sequence
                second_sequence = "-" + second_sequence
                gap_line = " " + gap_line
                i -= 1
            #Begin Gap
            elif top[i][j] == mid[i-1][j]-6:
                first_sequence = seq1[i-1] + first_sequence
                second_sequence = "-" + second_sequence
                gap_line = " " + gap_line
                i -= 1
                current_matrix = 1
            else:
                if i > 0:
                    i -= 1
                else:
                    j -= 1
        elif current_matrix == 1:
            #End of Top Gap
            if i != 0 and mid[i][j] == top[i][j]:
                current_matrix = 0
            #End of Bot Gap
            elif j != 0 and mid[i][j] == bot[i][j]:
                current_matrix = 2
            #Match
            elif mid[i][j] == mid[i-1][j-1] + 1:
                first_sequence = seq1[i-1] + first_sequence
                second_sequence = seq2[j-1] + second_sequence
                gap_line = "|" + gap_line
                j -= 1
                i -= 1
            #Mismatch
            elif mid[i][j] == mid[i-1][j-1] - 2:
                first_sequence = seq1[i-1] + first_sequence
                second_sequence = seq2[j-1] + second_sequence
                gap_line = " " + gap_line
                j -= 1
                i -= 1
        elif current_matrix == 2:
            #Extend Gap
            if bot[i][j] == bot[i][j-1] - 1:
                first_sequence = "-" + first_sequence
                second_sequence = seq2[j-1] + second_sequence
                gap_line = " " + gap_line
                j -= 1
            #Being Gap
            elif bot[i][j] == mid[i][j-1]-6:
                first_sequence = "-" + first_sequence
                second_sequence = seq2[j-1] + second_sequence
                gap_line = " " + gap_line
                j -= 1
                current_matrix = 1
        else:
            i -= 1
            j -= 1
    
    # print(current_matrix,i,j)
    print("\nThe alignment is:")
    for i in range(len(first_sequence)//70+1):
        if len(first_sequence) < (i+1)*70:
            print(first_sequence[i*70:].rstrip())
            print(gap_line[i*70:].rstrip())
            print(second_sequence[i*70:].rstrip())
        else:
            print(first_sequence[i*70:(i+1)*70].rstrip())
            print(gap_line[i*70:(i+1)*70].rstrip())
            print(second_sequence[i*70:(i+1)*70].rstrip())

def print_local_alignment(alignment_matrix,seq1,seq2):
    first_sequence = ""
    second_sequence = ""
    gap_line = ""
    max_i = 0
    max_j = 0
    current_max = 0
    for x in range(len(alignment_matrix)-1):
        for y in range(len(alignment_matrix[0])-1):
            if alignment_matrix[x][y] > current_max:
                max_i = x
                max_j = y
                current_max = alignment_matrix[x][y]
    
    i = max_i
    j = max_j
    while (i > 0 or j > 0) and alignment_matrix[i][j] != 0:
        # print("i:",i,"j:",j)
        # print(first_sequence)
        # print(gap_line)
        # print(second_sequence)
        current = alignment_matrix[i][j]
        if i > 0 and current - alignment_matrix[i-1][j] == -3:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = "-" + second_sequence
            gap_line = " " + gap_line
            i -= 1
            # print(1)
        elif current - alignment_matrix[i][j-1] == -3:
            first_sequence = "-" + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = " " + gap_line
            j -= 1
            # print(2)
        elif i > 0 and j > 0 and current - alignment_matrix[i-1][j-1] == 1:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = "|" + gap_line
            j -= 1
            i -= 1
            # print(3)
        elif i > 0 and j > 0 and current - alignment_matrix[i-1][j-1] == -2:
            first_sequence = seq1[i-1] + first_sequence
            second_sequence = seq2[j-1] + second_sequence
            gap_line = " " + gap_line
            j -= 1
            i -= 1
            # print(4)
        else:
            i -= 1
            j -= 1
            print(5)
    print("\nThe alignment is:")
    for i in range(len(first_sequence)//70+1):
        if len(first_sequence) < (i+1)*70:
            print(first_sequence[i*70:].rstrip())
            print(gap_line[i*70:].rstrip())
            print(second_sequence[i*70:].rstrip())
        else:
            print(first_sequence[i*70:(i+1)*70].rstrip())
            print(gap_line[i*70:(i+1)*70].rstrip())
            print(second_sequence[i*70:(i+1)*70].rstrip())

def print_alignment_matrix(matrix,fill=4):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            print(str(matrix[i][j]).zfill(fill), end =" ")
        print()