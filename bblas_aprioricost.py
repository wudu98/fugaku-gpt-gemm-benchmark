import sys
import csv
import os.path
import shutil

src_path = './bblas_src_aprioricost'

if not os.path.exists(src_path):
  os.mkdir(src_path)

# ask whether the argument is pointer
#  return
#   argument is pointer : True
#             not       : False
def pointer_(typ):
  if typ.find('*') >= 0:
    ret = True
  else:
    ret = False
  return ret

# ask whether the argument's category is local
#  return
#   argument is local : True
#             not     : False
def local_(cat):
  if cat.find('l') >= 0:
    ret = True
  else:
    ret = False
  return ret

# ask whether the argument's category is group
#  return
#   argument is group : True
#             not     : False
def group_(cat):
  if cat.find('g') >= 0:
    ret = True
  else:
    ret = False
  return ret

# ask whether fuction doesn't have return value
#  return
#   not have return value : True
#       have return value : Flase
def void_(r):
  if r.find('void') >= 0:
    ret = True
  else:
    ret = False
  return ret

# routine declaration
def routine_declaration(fw,row1,row2):
  fw.write('void ')
  fw.write(row1[0])
  c = 0
  for r in row2:
    if c == 0:
      fw.write('(')
      if '_batchf' in row1[0] :
        fw.write('const int group_size,')
      else :
        fw.write('const int group_count,const int *group_size,')
    # get category
    if c % 3 == 0:
      typ = r 
      cat = row2[c+2]
    # write word
    if c % 3 == 2:
      pass
    else:
      fw.write(r)
    
    if c % 3 == 0:
      if (pointer_(typ) and local_(cat)) or \
            ((not pointer_(typ)) and (group_(cat) or local_(cat))):
        fw.write('*')
      fw.write(' ')
    elif c % 3 == 1:
      fw.write(',')
      
    c = c + 1

  # additional argument for returning result
  if not void_(row1[1]):
    fw.write(row1[1])
    fw.write('*')
    fw.write(' ret_')
    fw.write(row1[2])
    fw.write(',')
    
  fw.write('int *info)')

  return

# routine call
def routine_call(fw,row1,row2):

  if not void_(row1[1]):
    fw.write('ret_')
    fw.write(row1[2])
    fw.write('[group_head[group_no]+local_no]')
    fw.write(' = ')
  
  fw.write(row1[2])
  
  c = 0
  for r in row2:
    if c % 3 == 0:
      if c == 0:
        fw.write('(')
        
      cat = row2[c+2]
      # for complex number
      if cat.find('g') >= 0:
        if cat.find('z') >= 0:
          fw.write('(double *)')
        elif cat.find('c') >= 0:
          fw.write('(float *)')

      # write name
      fw.write(row2[c+1])
      
      if cat.find('l') >= 0:
        fw.write('[group_head[group_no]+local_no]')
      elif cat.find('g') >= 0:
        # for complex number
        if cat.find('z') >= 0 or cat.find('c') >= 0:
          fw.write('+2*group_no')
        else:
          fw.write('[group_no]')
      if c == len(row2) - 3:
        fw.write(');')
      else:
        fw.write(',')
       
    c = c + 1
    
  fw.write('\n')

  return

def routine_fixed(fw_r,row1,row2):
  fw_r.write('  int group_count=1;\n')
  fw_r.write('  '+row1[0].replace('_batchf','_batch'))
  c = 0
  a ='(group_count,&group_size,'
  for r in row2:
    c = c + 1
    if(c == 2):
      a=a + r + ','
    if(c == 3):
      c = 0   

  # additional argument for returning result
  if not void_(row1[1]):
    a = a + 'ret_' + row1[2] +','

  a = a + 'info); \n'
  a=a.replace(',incx,',',&incx,')
  a=a.replace(',incy,',',&incy,')

  a=a.replace(',n,',',&n,')
  a=a.replace(',m,',',&m,')
  a=a.replace(',k,',',&k,')
  a=a.replace(',a1,',',&a1,')
  a=a.replace(',sb,',',&sb,')
  a=a.replace(',c1,',',&c1,')
  a=a.replace(',s1,',',&s1,')
  a=a.replace(',transa,',',&transa,')
  a=a.replace(',transb,',',&transb,')
  a=a.replace(',trans,',',&trans,')
  a=a.replace(',diag,',',&diag,')
  a=a.replace(',kl,',',&kl,')
  a=a.replace(',ku,',',&ku,')
  a=a.replace(',side,',',&side,')

  a=a.replace(',alpha,',',&alpha,')
  a=a.replace(',beta,',',&beta,')
  a=a.replace(',lda,',',&lda,')
  a=a.replace(',ldb,',',&ldb,')
  a=a.replace(',ldc,',',&ldc,')
  a=a.replace(',uplo,',',&uplo,')

  fw_r.write(a)
  fw_r.write('}\n')

def routine_my_sgemm(fw_r):

  UNROLL_K = 4
  NR_MAIN = 4

  SIMD_LANE = 16
  assert (SIMD_LANE == 16)

  CONST_UNROLL_LANE = 1
  assert (CONST_UNROLL_LANE == 1)

  DATA_BYTE = 4
  assert (DATA_BYTE == 4)

  RESERVED_REG_NUM = 6

  # print('M=%d, N=%d, K=%d' % (M, N, K))

  def micro_kernel_loop_asm(LOOP_ID, LAST_K_ID, LINES, COLS, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG):
      code_str = ""

      UNROLL_NR = 2
      if COLS % 2 != 0 :
        UNROLL_NR = 1
      if LOOP_ID == LAST_K_ID and WITH_BIAS_FLAG :
        UNROLL_NR = COLS

      A_odd_flag = (LOOP_ID // CONST_UNROLL_LANE) % 2
      B_odd_flag = ((LOOP_ID * COLS + VEC_REG_B_LEN) // COLS) % 2
      ptr_B_POS = (LOOP_ID * COLS + VEC_REG_B_LEN) % COLS 
      mod_simd_lane_loop_id = LOOP_ID % CONST_UNROLL_LANE

      vector_scroll_A = [[], []]
      vector_scroll_A[0] = [ vector_id_array_A[i] for i in range(LINES)]
      vector_scroll_A[1] = [ vector_id_array_A[(i+real_lines)%VEC_REG_A_LEN ] for i in range(real_lines) ]
      
      vector_scroll_B = [(i + LOOP_ID * COLS) % VEC_REG_B_LEN for i in range(COLS)]

      # Initializes the ABC Block pointer
      if REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID:
        if REG_BLOCK_TRANS_FLAG == 2 :
          code_str += f"    \"mov     x21, %[A]                 \\n\"\n"
          code_str += f"    \"add     %[B], %[B], #{real_cols * DATA_BYTE}                 \\n\"\n"
          code_str += f"    \"add     %[C], %[C], #{real_cols * DATA_BYTE}                 \\n\"\n"
          code_str += f"    \"mov     x24, %[C]                 \\n\"\n"
        code_str += f"    \"mov     x22, %[B]                   \\n\"\n"
        code_str += f"    \"add     x23, %[B], %[ldb], lsr #1               \\n\"\n" 
        code_str += f"    \"prfm    PLDL1KEEP, [x22, #64]              \\n\"\n"
        code_str += f"    \"prfm    PLDL1KEEP, [x23, #64]              \\n\"\n"
        B_odd_flag = 0

      for i in range(LINES*COLS//UNROLL_NR):
        line = i % LINES
        col = i // LINES
        
        # Main computing
        if FMA_CALCULATE_FLAG :
          if(LOOP_ID == 0 and LOOP_K_BEGIN_FLAG and (not WITH_BIAS_FLAG)):
            for j in range(UNROLL_NR):
              if(line < real_lines and SIMD_LANE*UNROLL_NR*col + SIMD_LANE*j < real_cols):
                if CONST_UNROLL_LANE == 1 :
                  code_str += f"    \"fmul    z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, z{vector_scroll_A[A_odd_flag][line]}.s             \\n\"\n"
                else :
                  code_str += f"    \"fmul    z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, z{vector_scroll_A[A_odd_flag][line]}.s[{mod_simd_lane_loop_id}]             \\n\"\n"
                
          else:
            for j in range(UNROLL_NR):
              if(line < real_lines and SIMD_LANE*UNROLL_NR*col + SIMD_LANE*j < real_cols):
                if A_odd_flag == 1 and ((LOOP_ID == LAST_K_ID and not WITH_BIAS_FLAG) or (not LOOP_ID == LAST_K_ID and mod_simd_lane_loop_id == CONST_UNROLL_LANE - 1)) :
                  ori_line = line
                  line = (line + VEC_REG_A_LEN % real_lines) % real_lines
                if CONST_UNROLL_LANE == 1 :
                  if(SIMD_LANE*UNROLL_NR*col + SIMD_LANE*(j+1) <= real_cols):
                    code_str += f"    \"fmla    z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p0/m, z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, z{vector_scroll_A[A_odd_flag][line]}.s             \\n\"\n"
                  else :
                    code_str += f"    \"fmla    z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p1/m, z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, z{vector_scroll_A[A_odd_flag][line]}.s             \\n\"\n"
                else : 
                  code_str += f"    \"fmla    z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, z{vector_scroll_A[A_odd_flag][line]}.s[{mod_simd_lane_loop_id}]             \\n\"\n"
                if A_odd_flag == 1 and ((LOOP_ID == LAST_K_ID and not WITH_BIAS_FLAG) or (not LOOP_ID == LAST_K_ID and mod_simd_lane_loop_id == CONST_UNROLL_LANE - 1)) :
                  line = ori_line
                
                
        # Store C
        if(STORE_C_FLAG and LOOP_ID == LAST_K_ID):
          for j in range(UNROLL_NR):
            tmp_col = UNROLL_NR*col + j
            if line < real_lines :
              if A_odd_flag == 1 and (not WITH_BIAS_FLAG) :
                ori_line = line
                line = (line + VEC_REG_A_LEN % real_lines) % real_lines
              if(SIMD_LANE*UNROLL_NR*col + SIMD_LANE*(j+1) <= real_cols):
                if tmp_col == 0 :
                  code_str += f"    \"st1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p0, [x{RESERVED_REG_NUM+line}]           \\n\"\n"
                else :
                  code_str += f"    \"st1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p0, [x{RESERVED_REG_NUM+line}, #{tmp_col}, mul vl]           \\n\"\n"
              else:
                if tmp_col == 0 :
                  code_str += f"    \"st1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p1, [x{RESERVED_REG_NUM+line}]           \\n\"\n"
                else :
                  code_str += f"    \"st1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p1, [x{RESERVED_REG_NUM+line}, #{tmp_col}, mul vl]           \\n\"\n"
              if A_odd_flag == 1 and (not WITH_BIAS_FLAG) :
                line = ori_line
        
        if LOOP_K_END_FLAG and LOOP_ID == LAST_K_ID:
          continue

        if not WITH_BIAS_FLAG:
          # Get next block C address
          if (REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID and line == LINES - 1 and col == COLS//UNROLL_NR - 1):
            for j in range(next_lines):
              if (j == 0):
                code_str += f"    \"mov     x{RESERVED_REG_NUM}, x24    \\n\"\n"
              elif(j == 1):
                code_str += f"    \"add     x{RESERVED_REG_NUM+1}, x24, %[ldc]     \\n\"\n"
              else:
                code_str += f"    \"add     x{RESERVED_REG_NUM+j}, x{RESERVED_REG_NUM+j-2}, %[ldc], lsl #1    \\n\"\n"
        else:
          # Get next block C address
          if (REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID):
            if line < next_lines :
              if (line == 0):
                code_str += f"    \"mov     x{RESERVED_REG_NUM}, x24    \\n\"\n"
              elif(line == 1):
                code_str += f"    \"add     x{RESERVED_REG_NUM+1}, x24, %[ldc]     \\n\"\n"
              else:
                code_str += f"    \"add     x{RESERVED_REG_NUM+line}, x{RESERVED_REG_NUM+line-2}, %[ldc], lsl #1    \\n\"\n"
          # Load next block C in vector register
          if REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID:
            for j in range(UNROLL_NR):
              tmp_col = UNROLL_NR*col + j
              if(line < next_lines and SIMD_LANE*UNROLL_NR*col + SIMD_LANE*j < next_cols):
                if tmp_col == 0 :
                  code_str += f"    \"ld1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p0/z, [x{RESERVED_REG_NUM+line}]           \\n\"\n"
                else :
                  code_str += f"    \"ld1w     z{VEC_REG_A_LEN + VEC_REG_B_LEN + line*COLS + col*UNROLL_NR + j}.s, p0/z, [x{RESERVED_REG_NUM+line}, #{tmp_col}, mul vl]           \\n\"\n"

        # Get next block A address
        if (REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID and line == 0 and col == 0):
          for j in range(next_lines):
            if (j == 0):
              code_str += f"    \"mov     x{RESERVED_REG_NUM+LINES}, x21    \\n\"\n"
            elif(j == 1):
              code_str += f"    \"add     x{RESERVED_REG_NUM+LINES+1}, x21, %[lda]    \\n\"\n"
            else:
              code_str += f"    \"add     x{RESERVED_REG_NUM+LINES+j}, x{RESERVED_REG_NUM+LINES+j-2}, %[lda], lsl #1    \\n\"\n"

        # Load next A in vector register
        if not REG_BLOCK_TRANS_FLAG :
          # Sequence load next A
          # The code contains the permutation operation, registers used to scroll A to improve performance
          # Corresponding to the Main computing
          if ((LAST_K_ID == -1 or LOOP_ID < (LAST_K_ID - LAST_K_ID%4)) and line == 0 and col == 0):
            ori_line = line
            for line in range(real_lines):
              if((CONST_UNROLL_LANE - 1 == 0 or mod_simd_lane_loop_id == line % (CONST_UNROLL_LANE - 1)) and (line >= real_lines - VEC_REG_A_LEN % real_lines or 2 * real_lines <= VEC_REG_A_LEN)):
                if A_odd_flag == 0:
                  line = (line + VEC_REG_A_LEN % real_lines) % real_lines
                if CONST_UNROLL_LANE == 1 :
                  code_str += f"    \"ld1rw      z{vector_scroll_A[A_odd_flag^1][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
                else :
                  code_str += f"    \"ld1rqw      z{vector_scroll_A[A_odd_flag^1][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
                code_str += f"    \"add      x{RESERVED_REG_NUM+LINES+line}, x{RESERVED_REG_NUM+LINES+line}, #{CONST_UNROLL_LANE * DATA_BYTE}     \\n\"\n"
            line = ori_line

          if((LAST_K_ID == -1 or LOOP_ID < (LAST_K_ID - LAST_K_ID%4)) and mod_simd_lane_loop_id == CONST_UNROLL_LANE - 1 and line < real_lines and col == (real_cols+SIMD_LANE-1)//SIMD_LANE//UNROLL_NR - 1):
            if (2 * real_lines > VEC_REG_A_LEN and line < real_lines - VEC_REG_A_LEN % real_lines):
              if A_odd_flag == 0:
                ori_line = line
                line = (line + VEC_REG_A_LEN % real_lines) % real_lines
              if CONST_UNROLL_LANE == 1 :
                code_str += f"    \"ld1rw     z{vector_scroll_A[A_odd_flag^1][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
              else :
                code_str += f"    \"ld1rqw     z{vector_scroll_A[A_odd_flag^1][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
              code_str += f"    \"add      x{RESERVED_REG_NUM+LINES+line}, x{RESERVED_REG_NUM+LINES+line}, #{CONST_UNROLL_LANE * DATA_BYTE}     \\n\"\n"
              if A_odd_flag == 0:
                line = ori_line
        else :
          if not WITH_BIAS_FLAG:
            # Load next block A
            if(LOOP_ID == LAST_K_ID and line < next_lines and col == (real_cols+SIMD_LANE-1)//SIMD_LANE//UNROLL_NR - 1) :
              if CONST_UNROLL_LANE == 1 :
                code_str += f"    \"ld1rw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
              else :
                code_str += f"    \"ld1rqw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
              code_str += f"    \"add      x{RESERVED_REG_NUM+LINES+line}, x{RESERVED_REG_NUM+LINES+line}, #{CONST_UNROLL_LANE * DATA_BYTE}     \\n\"\n"
          else:
            # Load next block A
            if(LOOP_ID == LAST_K_ID and line < next_lines and col == (real_cols+SIMD_LANE-1)//SIMD_LANE//UNROLL_NR - 1) :
              if A_odd_flag == 0 or line >= real_lines - VEC_REG_A_LEN % real_lines:
                if CONST_UNROLL_LANE == 1 :
                  code_str += f"    \"ld1rw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
                else :
                  code_str += f"    \"ld1rqw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
                code_str += f"    \"add      x{RESERVED_REG_NUM+LINES+line}, x{RESERVED_REG_NUM+LINES+line}, #{CONST_UNROLL_LANE * DATA_BYTE}     \\n\"\n"

        # Sequence Load next B in vector register
        if (line == LINES - 1):
          for j in range(UNROLL_NR):
            if(((not LOOP_ID == LAST_K_ID) or (LOOP_ID == LAST_K_ID and COLS == VEC_REG_B_LEN)) and SIMD_LANE*UNROLL_NR*col + SIMD_LANE*j < real_cols):
              if (LOOP_ID == LAST_K_ID - 1 and UNROLL_NR*col + j >= 2 * COLS - VEC_REG_B_LEN):
                continue
              if ptr_B_POS == 0 :
                code_str += f"    \"ld1w     z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, p0/z, [x{register_scroll_B[B_odd_flag]}]             \\n\"\n"
              else :
                code_str += f"    \"ld1w     z{vector_id_array_B[vector_scroll_B[col*UNROLL_NR + j]]}.s, p0/z, [x{register_scroll_B[B_odd_flag]}, #{ptr_B_POS}, mul vl]             \\n\"\n"
              # Get next B address
              if ptr_B_POS == COLS - 1:
                ptr_B_POS = 0
                code_str += f"    \"add     x{register_scroll_B[B_odd_flag]}, x{register_scroll_B[B_odd_flag]}, %[ldb]              \\n\"\n"
                code_str += f"    \"prfm    PLDL1KEEP, [x{register_scroll_B[B_odd_flag]}, #64]              \\n\"\n"
                B_odd_flag ^= 1
              else:
                ptr_B_POS += 1 
      
      # Extra operations ensure that load next block A works correctly
      if REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID and WITH_BIAS_FLAG:
        for line in range(next_lines):
          if not (A_odd_flag == 0 or line >= real_lines - VEC_REG_A_LEN % real_lines):
            if CONST_UNROLL_LANE == 1 :
              code_str += f"    \"ld1rw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
            else :
              code_str += f"    \"ld1rqw     z{vector_scroll_A[0][line]}.s, p0/z, [x{RESERVED_REG_NUM+LINES+line}]    \\n\"\n"
            code_str += f"    \"add      x{RESERVED_REG_NUM+LINES+line}, x{RESERVED_REG_NUM+LINES+line}, #{CONST_UNROLL_LANE * DATA_BYTE}     \\n\"\n"
      
      # Extra operations ensure that Load next block B works correctly
      if REG_BLOCK_TRANS_FLAG and LOOP_ID == LAST_K_ID and (not COLS == VEC_REG_B_LEN):
        vector_scroll_B = [i for i in range(VEC_REG_B_LEN)]
        ptr_B_POS = 0
        for j in range(VEC_REG_B_LEN):
          if ptr_B_POS == 0 :
            code_str += f"    \"ld1w     z{vector_id_array_B[vector_scroll_B[j]]}.s, p0/z, [x{register_scroll_B[B_odd_flag]}]             \\n\"\n"
          else :
            code_str += f"    \"ld1w     z{vector_id_array_B[vector_scroll_B[j]]}.s, p0/z, [x{register_scroll_B[B_odd_flag]}, #{ptr_B_POS}, mul vl]             \\n\"\n"
          if ptr_B_POS == COLS - 1:
            ptr_B_POS = 0
            code_str += f"    \"add     x{register_scroll_B[B_odd_flag]}, x{register_scroll_B[B_odd_flag]}, %[ldb]              \\n\"\n"
            code_str += f"    \"prfm    PLDL1KEEP, [x{register_scroll_B[B_odd_flag]}, #64]              \\n\"\n"
            B_odd_flag ^= 1
          else:
            ptr_B_POS += 1 

      return code_str

  def UNROLL_LOOP_ID(K, UNROLL_K):
    BEGIN_LOOP = 1
    assert (BEGIN_LOOP == 1)
    EDGE_BEGIN_LOOP = 1
    if K % UNROLL_K > CONST_UNROLL_LANE :
      EDGE_BEGIN_LOOP = (K % UNROLL_K) - (K % CONST_UNROLL_LANE)
    elif K % UNROLL_K == 0 and UNROLL_K > CONST_UNROLL_LANE :
      EDGE_BEGIN_LOOP = UNROLL_K - CONST_UNROLL_LANE
    return BEGIN_LOOP, EDGE_BEGIN_LOOP

  def compile_time_for_init_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 1
      FMA_CALCULATE_FLAG = 0
      STORE_C_FLAG = 0
      WITH_BIAS_FLAG = with_bias

      code_str = ""
      code_str += micro_kernel_loop_asm(-1, -1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str

  def compile_time_for_loop_k_begin_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias):
      LOOP_K_BEGIN_FLAG = 1
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 0
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 0
      WITH_BIAS_FLAG = with_bias

      MAIN_K_LOOP_BEGIN = 0
      MAIN_K_LOOP_END, _ = UNROLL_LOOP_ID(K, UNROLL_K)

      code_str = f""

      tmp_LINES = real_lines
      cnt = 0
      while(tmp_LINES != 0):
        if (tmp_LINES % 2 != 0):
          if cnt == 0:
            code_str += f"    \"add     x21, x21, %[lda]               \\n\"\n"
            code_str += f"    \"add     x24, x24, %[ldc]               \\n\"\n"
          else:
            code_str += f"    \"add     x21, x21, %[lda], lsl #{cnt}               \\n\"\n"
            code_str += f"    \"add     x24, x24, %[ldc], lsl #{cnt}               \\n\"\n"
        tmp_LINES = tmp_LINES // 2
        cnt += 1

      for LOOP_ID in range(MAIN_K_LOOP_BEGIN, MAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, -1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str

  def compile_time_for_loop_k_main_body_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 0
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 0
      WITH_BIAS_FLAG = 0

      MAIN_K_LOOP_BEGIN, _ = UNROLL_LOOP_ID(K, UNROLL_K)
      MAIN_K_LOOP_END = UNROLL_K + MAIN_K_LOOP_BEGIN

      code_str = f""

      for LOOP_ID in range(MAIN_K_LOOP_BEGIN, MAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, -1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str


  def compile_time_for_loop_k_remain_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 0
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 0
      WITH_BIAS_FLAG = 0

      REMAIN_K_LOOP_BEGIN, REMAIN_K_LOOP_END = UNROLL_LOOP_ID(K, UNROLL_K)

      code_str = f""

      for LOOP_ID in range(REMAIN_K_LOOP_BEGIN, REMAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, -1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str
      
      
  def compile_time_for_m_dim_micro_kernel_pipeline_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 1
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 1
      WITH_BIAS_FLAG = with_bias
      
      _, REMAIN_K_LOOP_BEGIN = UNROLL_LOOP_ID(K, UNROLL_K)
      REMAIN_K_LOOP_END = UNROLL_K if K % UNROLL_K == 0 else K % UNROLL_K
      
      code_str = f""

      for line in range(real_lines):
        code_str += f"    \"prfm    PSTL1KEEP, [x{RESERVED_REG_NUM+line}, #64]              \\n\"\n"

      for LOOP_ID in range(REMAIN_K_LOOP_BEGIN, REMAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, REMAIN_K_LOOP_END-1, LINES, COLS, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      # When K module UNROLL_K remainder 1, no calculation, direct store
      if (REMAIN_K_LOOP_BEGIN == REMAIN_K_LOOP_END):
        FMA_CALCULATE_FLAG = 0
        code_str += micro_kernel_loop_asm(-1, -1, LINES, COLS, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str

  def compile_time_for_n_dim_micro_kernel_pipeline_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 0
      REG_BLOCK_TRANS_FLAG = 2
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 1
      WITH_BIAS_FLAG = with_bias

      _, REMAIN_K_LOOP_BEGIN = UNROLL_LOOP_ID(K, UNROLL_K)
      REMAIN_K_LOOP_END = UNROLL_K if K % UNROLL_K == 0 else K % UNROLL_K
      
      code_str = f""

      for line in range(real_lines):
        code_str += f"    \"prfm    PSTL1KEEP, [x{RESERVED_REG_NUM+line}, #64]              \\n\"\n"

      for LOOP_ID in range(REMAIN_K_LOOP_BEGIN, REMAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, REMAIN_K_LOOP_END-1, LINES, COLS, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      # When K module UNROLL_K remainder 1, no calculation, direct store
      if (REMAIN_K_LOOP_BEGIN == REMAIN_K_LOOP_END):
        FMA_CALCULATE_FLAG = 0
        code_str += micro_kernel_loop_asm(-1, -1, LINES, COLS, real_lines, real_cols, next_lines, next_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str

  def compile_time_for_loop_k_end_func_asm(LINES, COLS, K, UNROLL_K, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B):
      LOOP_K_BEGIN_FLAG = 0
      LOOP_K_END_FLAG = 1
      REG_BLOCK_TRANS_FLAG = 0
      FMA_CALCULATE_FLAG = 1
      STORE_C_FLAG = 1
      WITH_BIAS_FLAG = 0

      _, REMAIN_K_LOOP_BEGIN = UNROLL_LOOP_ID(K, UNROLL_K)
      REMAIN_K_LOOP_END = UNROLL_K if K % UNROLL_K == 0 else K % UNROLL_K
      
      code_str = f""

      for line in range(real_lines):
        code_str += f"    \"prfm    PSTL1KEEP, [x{RESERVED_REG_NUM+line}, #64]              \\n\"\n"

      for LOOP_ID in range(REMAIN_K_LOOP_BEGIN, REMAIN_K_LOOP_END):
        code_str += micro_kernel_loop_asm(LOOP_ID, REMAIN_K_LOOP_END-1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      # When K module UNROLL_K remainder 1, no calculation, direct store
      if (REMAIN_K_LOOP_BEGIN == REMAIN_K_LOOP_END):
        FMA_CALCULATE_FLAG = 0
        code_str += micro_kernel_loop_asm(-1, -1, LINES, COLS, real_lines, real_cols, real_lines, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, LOOP_K_BEGIN_FLAG, LOOP_K_END_FLAG, REG_BLOCK_TRANS_FLAG, FMA_CALCULATE_FLAG, STORE_C_FLAG, WITH_BIAS_FLAG)

      return code_str


  def m_dim_func_asm(MR_MAIN, MR_MAIN_LOOPS, MR_REMAIN, MR_REMAIN_LOOPS, NR, K, UNROLL_K, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias):
      
      Main_K_loop_flag = (K > UNROLL_K) # The K-dim main operation needs to loop through
      Main_K_loop_times = (K + UNROLL_K - 1) // UNROLL_K
      
      code_str = f""

      if MR_MAIN_LOOPS : # Enter the M-dim main operation
        if MR_MAIN_LOOPS > 1 : # Cyclic M-dim main operation
          code_str += f"    \"mov     x26, #{MR_MAIN_LOOPS}                   \\n\"\n"
          code_str += f"    \"b       1f                                 \\n\"\n"
          code_str += f"  \"2:                                 \\n\"\n"
          code_str += f"    \"subs    x26, x26, #1                            \\n\"\n"
          code_str += compile_time_for_loop_k_remain_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)
          code_str += f"    \"beq     3f              \\n\"\n"
          code_str += compile_time_for_m_dim_micro_kernel_pipeline_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          code_str += f"  \"1:                                 \\n\"\n"

        # K-dim main operation
        if Main_K_loop_flag : 
          code_str += f"    \"mov     x27, #{Main_K_loop_times}                   \\n\"\n"
          code_str += f"    \"subs    x27, x27, #1                            \\n\"\n"
          code_str += compile_time_for_loop_k_begin_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          code_str += f"    \"b       4f                                 \\n\"\n"
          code_str += f"  \"5:                                 \\n\"\n"
          code_str += compile_time_for_loop_k_main_body_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)
          code_str += f"  \"4:                                 \\n\"\n"
          if MR_MAIN_LOOPS > 1 :
            code_str += f"    \"beq     2b                       \\n\"\n"
          else:
            code_str += f"    \"beq     3f                       \\n\"\n"
          code_str += f"    \"subs    x27, x27, #1                            \\n\"\n"
          code_str += f"    \"b       5b                                 \\n\"\n"
        else:
          code_str += compile_time_for_loop_k_begin_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          if MR_MAIN_LOOPS > 1 :
            code_str += f"    \"b       2b                       \\n\"\n"

        if MR_MAIN_LOOPS > 1 or Main_K_loop_flag:
          code_str += f"  \"3:                                 \\n\"\n"

        if not MR_MAIN_LOOPS > 1 :
          code_str += compile_time_for_loop_k_remain_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)

      if MR_REMAIN_LOOPS : # Enter the M-dim remain operation
        if MR_MAIN_LOOPS : # Cyclic M-dim remain operation
          code_str += compile_time_for_m_dim_micro_kernel_pipeline_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_MAIN, real_cols, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
        if MR_REMAIN_LOOPS > 1 :
          code_str += f"    \"mov     x26, #{MR_REMAIN_LOOPS}                   \\n\"\n"
          code_str += f"    \"b       1f                                 \\n\"\n"
          code_str += f"  \"2:                                 \\n\"\n"
          code_str += f"    \"subs    x26, x26, #1                            \\n\"\n"
          code_str += compile_time_for_loop_k_remain_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)
          code_str += f"    \"beq     3f              \\n\"\n"
          code_str += compile_time_for_m_dim_micro_kernel_pipeline_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          code_str += f"  \"1:                                 \\n\"\n"

        # K-dim main operation
        if Main_K_loop_flag : 
          code_str += f"    \"mov     x27, #{Main_K_loop_times}                   \\n\"\n"
          code_str += f"    \"subs    x27, x27, #1                            \\n\"\n"
          code_str += compile_time_for_loop_k_begin_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          code_str += f"    \"b       4f                                 \\n\"\n"
          code_str += f"  \"5:                                 \\n\"\n"
          code_str += compile_time_for_loop_k_main_body_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)
          code_str += f"  \"4:                                 \\n\"\n"
          if MR_REMAIN_LOOPS > 1 :
            code_str += f"    \"beq     2b                       \\n\"\n"
          else:
            code_str += f"    \"beq     3f                       \\n\"\n"
          code_str += f"    \"subs    x27, x27, #1                            \\n\"\n"
          code_str += f"    \"b       5b                                 \\n\"\n"
        else:
          code_str += compile_time_for_loop_k_begin_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          if MR_REMAIN_LOOPS > 1 :
            code_str += f"    \"b       2b                       \\n\"\n"

        if MR_REMAIN_LOOPS > 1 or Main_K_loop_flag:
          code_str += f"  \"3:                                 \\n\"\n"

        if not MR_REMAIN_LOOPS > 1 :
          code_str += compile_time_for_loop_k_remain_func_asm(MR_MAIN, NR, K, UNROLL_K, MR_REMAIN, real_cols, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)
          
      return code_str

  def n_dim_func_asm(REMAIN_N, K, UNROLL_K, NR, NR_LOOPS, MR_MAIN, MR_MAIN_LOOPS, MR_REMAIN, MR_REMAIN_LOOPS, with_bias) :

      VEC_REG_B_LEN = NR if K <= 16 else max(4, NR)
      if NR == 6 :
        VEC_REG_B_LEN = NR if K <= 32 else 8

      VEC_REG_A_LEN = MR_MAIN if K <= 16 else min(32 - MR_MAIN*NR - VEC_REG_B_LEN, min(2*MR_MAIN, 8))

      vector_id_array_A = []
      for i in range(0, VEC_REG_A_LEN):
        vector_id_array_A.append(i)
      
      vector_id_array_B = []
      for i in range(VEC_REG_A_LEN, VEC_REG_A_LEN+VEC_REG_B_LEN):
        vector_id_array_B.append(i)

      register_scroll_B = [22, 23]

      Edge_N_flag = SIMD_LANE * NR * NR_LOOPS > REMAIN_N
      Edge_N = REMAIN_N % (SIMD_LANE * NR)
      if Edge_N_flag : 
        NR_LOOPS -= 1
      Main_N_flag = 0 if NR_LOOPS == 0 else 1
      
      lines_branch_1 = MR_MAIN if MR_MAIN_LOOPS else MR_REMAIN
      lines_branch_2 = MR_MAIN if not MR_REMAIN_LOOPS else MR_REMAIN
      cols_branch_1 = SIMD_LANE*NR if Main_N_flag else Edge_N
      cols_branch_2 = SIMD_LANE*NR if not Edge_N_flag else Edge_N

      code_str = f""
      code_str += f"    \"ptrue     p0.s                  \\n\"\n"
      code_str += f"    \"mov       x28, #{SIMD_LANE if REMAIN_N % SIMD_LANE == 0 else REMAIN_N % SIMD_LANE}                  \\n\"\n"
      code_str += f"    \"whilelt   p1.s, xzr, x28                  \\n\"\n"
      code_str += compile_time_for_init_func_asm(MR_MAIN, NR, K, UNROLL_K, lines_branch_1, cols_branch_1, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)

      if Main_N_flag : # Enter the N-dim main operation
        if NR_LOOPS > 1 : # Cyclic N-dim main operation
          code_str += f"    \"mov     x25, #{NR_LOOPS}                   \\n\"\n"
          code_str += f"    \"b       6f                                 \\n\"\n"
          code_str += f"  \"0:                                 \\n\"\n"
          code_str += f"    \"subs    x25, x25, #1                            \\n\"\n"
          code_str += f"    \"beq     7f                       \\n\"\n" 
          code_str += compile_time_for_n_dim_micro_kernel_pipeline_func_asm(MR_MAIN, NR, K, UNROLL_K, lines_branch_2, SIMD_LANE*NR, lines_branch_1, SIMD_LANE*NR, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
          code_str += f"  \"6:                                 \\n\"\n"
        
        code_str += m_dim_func_asm(MR_MAIN, MR_MAIN_LOOPS, MR_REMAIN, MR_REMAIN_LOOPS, NR, K, UNROLL_K, SIMD_LANE*NR, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)

        if NR_LOOPS > 1 : 
          code_str += f"    \"b       0b                                 \\n\"\n"
          code_str += f"  \"7:                                 \\n\"\n"

      if Edge_N_flag : # Enter the N-dim remain operation
        if Main_N_flag : # Cyclic N-dim remain operation
          code_str += compile_time_for_n_dim_micro_kernel_pipeline_func_asm(MR_MAIN, NR, K, UNROLL_K, lines_branch_2, SIMD_LANE*NR, lines_branch_1, Edge_N, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)
        code_str += m_dim_func_asm(MR_MAIN, MR_MAIN_LOOPS, MR_REMAIN, MR_REMAIN_LOOPS, NR, K, UNROLL_K, Edge_N, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B, with_bias)

      code_str += compile_time_for_loop_k_end_func_asm(MR_MAIN, NR, K, UNROLL_K, lines_branch_2, cols_branch_2, vector_id_array_A, VEC_REG_A_LEN, vector_id_array_B, VEC_REG_B_LEN, register_scroll_B)

      return code_str

  def NRSA(N, NR_MAIN):
      CEIL_NC = (N + SIMD_LANE - 1) // SIMD_LANE
      NR_REMAIN = CEIL_NC % NR_MAIN
      NR_MAIN_LOOPS = CEIL_NC // NR_MAIN
      NR_REMAIN_LOOPS = 1 if NR_REMAIN else 0
      if NR_MAIN == 3 :
        if NR_REMAIN == 1 and NR_MAIN_LOOPS >= 1 :
          NR_MAIN_LOOPS -= 1
          NR_REMAIN = 4
        elif NR_REMAIN == 2 and NR_MAIN_LOOPS >= 1 :
          NR_MAIN_LOOPS -= 1
          NR_REMAIN = 5
      elif NR_MAIN == 4 :
        if NR_REMAIN == 1 and NR_MAIN_LOOPS >= 1 :
          NR_MAIN_LOOPS -= 1
          NR_REMAIN = 5
        elif NR_REMAIN == 2 :
          if NR_MAIN_LOOPS >= 2 :
            NR_MAIN_LOOPS -= 2
            NR_REMAIN = 5
            NR_REMAIN_LOOPS = 2
          elif NR_MAIN_LOOPS >= 1 :
            NR_MAIN_LOOPS -= 1
            NR_REMAIN = 6
        elif NR_REMAIN == 3 and NR_MAIN_LOOPS >= 3 :
          NR_MAIN_LOOPS -= 3
          NR_REMAIN = 5
          NR_REMAIN_LOOPS = 3
      elif NR_MAIN == 5 :
        if NR_REMAIN == 1 :
          if NR_MAIN_LOOPS >= 3 :
            NR_MAIN_LOOPS -= 3
            NR_REMAIN = 4
            NR_REMAIN_LOOPS = 4
          elif NR_MAIN_LOOPS >= 1 :
            NR_MAIN_LOOPS -= 1
            NR_REMAIN = 6
        elif NR_REMAIN == 2 and NR_MAIN_LOOPS >= 2 :
          NR_MAIN_LOOPS -= 2
          NR_REMAIN = 4
          NR_REMAIN_LOOPS = 3
        elif NR_REMAIN == 3 and NR_MAIN_LOOPS >= 1 :
          NR_MAIN_LOOPS -= 1
          NR_REMAIN = 4
          NR_REMAIN_LOOPS = 2

      return NR_MAIN_LOOPS, NR_REMAIN, NR_REMAIN_LOOPS

  def MRSA(M, NR):
      MR_MAIN = min(6, (32 - max(4, NR)) // (NR + 1))
      MR_REMAIN = M % MR_MAIN
      MR_MAIN_LOOPS = M // MR_MAIN
      MR_REMAIN_LOOPS = 1 if MR_REMAIN else 0
      if MR_MAIN == 5 :
        if MR_REMAIN == 1 :
          if MR_MAIN_LOOPS >= 3 :
            MR_MAIN_LOOPS -= 3
            MR_REMAIN = 4
            MR_REMAIN_LOOPS = 4
          elif MR_MAIN_LOOPS >= 1 :
            MR_MAIN_LOOPS -= 1
            MR_REMAIN = 3
            MR_REMAIN_LOOPS = 2
        elif MR_REMAIN == 2 and MR_MAIN_LOOPS >= 2 :
          MR_MAIN_LOOPS -= 2
          MR_REMAIN = 4
          MR_REMAIN_LOOPS = 3
        elif MR_REMAIN == 3 and MR_MAIN_LOOPS >= 1 :
          MR_MAIN_LOOPS -= 1
          MR_REMAIN = 4
          MR_REMAIN_LOOPS = 2
      elif MR_MAIN == 4 :
        if MR_REMAIN == 1 and MR_MAIN_LOOPS >= 2 :
          MR_MAIN_LOOPS -= 2
          MR_REMAIN = 3
          MR_REMAIN_LOOPS = 3
        elif MR_REMAIN == 2 and MR_MAIN_LOOPS >= 1 :
            MR_MAIN_LOOPS -= 1
            MR_REMAIN = 3
            MR_REMAIN_LOOPS = 2
      elif MR_MAIN == 3 and MR_REMAIN == 1 and MR_MAIN_LOOPS >= 1 :
        MR_MAIN_LOOPS -= 1
        MR_REMAIN = 2
        MR_REMAIN_LOOPS = 2
      
      return MR_MAIN, MR_MAIN_LOOPS, MR_REMAIN, MR_REMAIN_LOOPS

  def RBSA(M, N, NR_MAIN):
      NR_MAIN_LOOPS, NR_REMAIN, NR_REMAIN_LOOPS = NRSA(N, NR_MAIN)
      NR_MAIN_MR_MAIN, NR_MAIN_MR_MAIN_LOOPS, NR_MAIN_MR_REMAIN, NR_MAIN_MR_REMAIN_LOOPS = MRSA(M, NR_MAIN) if NR_MAIN_LOOPS else (0,0,0,0)
      NR_REMAIN_MR_MAIN, NR_REMAIN_MR_MAIN_LOOPS, NR_REMAIN_MR_REMAIN, NR_REMAIN_MR_REMAIN_LOOPS = MRSA(M, NR_REMAIN) if NR_REMAIN_LOOPS else (0,0,0,0)

      return NR_MAIN_LOOPS, NR_REMAIN, NR_REMAIN_LOOPS, NR_MAIN_MR_MAIN, NR_MAIN_MR_MAIN_LOOPS, NR_MAIN_MR_REMAIN, NR_MAIN_MR_REMAIN_LOOPS, NR_REMAIN_MR_MAIN, NR_REMAIN_MR_MAIN_LOOPS, NR_REMAIN_MR_REMAIN, NR_REMAIN_MR_REMAIN_LOOPS

  def laf_asm_code(M, N, K, UNROLL_K = 8, NR_MAIN = 4, with_bias = 0):

      assert (UNROLL_K % (2*CONST_UNROLL_LANE) == 0)
      assert (UNROLL_K >= 4)
      assert (NR_MAIN == 3 or NR_MAIN == 4 or NR_MAIN == 5)

      NR_MAIN_LOOPS, NR_REMAIN, NR_REMAIN_LOOPS, NR_MAIN_MR_MAIN, NR_MAIN_MR_MAIN_LOOPS, NR_MAIN_MR_REMAIN, NR_MAIN_MR_REMAIN_LOOPS, NR_REMAIN_MR_MAIN, NR_REMAIN_MR_MAIN_LOOPS, NR_REMAIN_MR_REMAIN, NR_REMAIN_MR_REMAIN_LOOPS = RBSA(M, N, NR_MAIN)

      code_str = ""
      code_str += f"""
    asm volatile(
      "prfm    PLDL1KEEP, [%[A], #64]              \\n"
      "prfm    PLDL1KEEP, [%[B], #64]              \\n"
      "lsl	   %[lda], %[lda], #2              \\n"
      "lsl	   %[ldb], %[ldb], #3              \\n"
      "lsl	   %[ldc], %[ldc], #2              \\n"
      "mov     x21, %[A]                 \\n"
      "mov     x24, %[C]                 \\n"
  """

      if NR_MAIN_LOOPS :
        code_str += n_dim_func_asm(min(N, SIMD_LANE * NR_MAIN * NR_MAIN_LOOPS), K, UNROLL_K, NR_MAIN, NR_MAIN_LOOPS, NR_MAIN_MR_MAIN, NR_MAIN_MR_MAIN_LOOPS, NR_MAIN_MR_REMAIN, NR_MAIN_MR_REMAIN_LOOPS, with_bias)

      if NR_REMAIN_LOOPS :
        if NR_MAIN_LOOPS:
          code_str += f"    \"mov     x21, %[A]                 \\n\"\n"
          code_str += f"    \"add     %[B], %[B], #{NR_MAIN * SIMD_LANE * 4}                 \\n\"\n"
          code_str += f"    \"add     %[C], %[C], #{NR_MAIN * SIMD_LANE * 4}                 \\n\"\n"
          code_str += f"    \"mov     x24, %[C]                 \\n\"\n"
        code_str += n_dim_func_asm(N - SIMD_LANE * NR_MAIN * NR_MAIN_LOOPS, K, UNROLL_K, NR_REMAIN, NR_REMAIN_LOOPS, NR_REMAIN_MR_MAIN, NR_REMAIN_MR_MAIN_LOOPS, NR_REMAIN_MR_REMAIN, NR_REMAIN_MR_REMAIN_LOOPS, with_bias)

      code_str += f"""
      : [A]"=r"(A),
        [B]"=r"(B),
        [C]"=r"(C),
        [lda]"=r"(lda),
        [ldb]"=r"(ldb),
        [ldc]"=r"(ldc)
      : "0"(A),
        "1"(B),
        "2"(C),
        "3"(lda),
        "4"(ldb),
        "5"(ldc)
      : "cc", "memory" """
      for i in range(6, 29):
        if i == 18 :
          continue
        else :
          code_str += f", \"x{i}\""
      code_str += f"\n                      "
      for i in range(32):
        code_str += f", \"z{i}\""
      code_str +=  f"""
    );
  """
      return code_str

  def xsmm_asm_armv8_code(UNROLL_K, NR_MAIN):
      MC = 64
      NC = 64
      KC = 128

      """Emit C code for gemm impl."""
      cc_code = f"""
  #ifndef __MY_SGEMM_KERNEL_H
  #define __MY_SGEMM_KERNEL_H
  #endif

void small_gemm_fixmn(const float *A, const float *B, float *C, int lda, int ldb, int ldc) {{
"""
      cc_code += laf_asm_code(MC, NC, KC, UNROLL_K, NR_MAIN, with_bias = 0)
      cc_code += f"""
}}

void small_gemm_fixmn_with_bias(const float *A, const float *B, float *C, int lda, int ldb, int ldc) {{
"""
      cc_code += laf_asm_code(MC, NC, KC, UNROLL_K, NR_MAIN, with_bias = 1)
      cc_code += f"""
}}

void sl_cblas_my_sgemm(const CBLAS_LAYOUT lay, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float  *A,
                 const int lda, const float  *B, const int ldb,
                 const float beta, float  *C, const int ldc) {{
  for (int j=0; j<N; j+={NC})
    for (int k=0; k<K; k+={KC})
      if (k == 0) {{
        for (int i=0; i<M; i+={MC})
          small_gemm_fixmn(A + i * lda + k, B + k * ldb + j, C + i * ldc + j, lda, ldb, ldc);
      }}
      else{{
        for (int i=0; i<M; i+={MC})
          small_gemm_fixmn_with_bias(A + i * lda + k, B + k * ldb + j, C + i * ldc + j, lda, ldb, ldc);
      }}
  }}
  
  """
      return cc_code

  fw_r.write(xsmm_asm_armv8_code(UNROLL_K, NR_MAIN))

  return

def routine_part1(fw_r):
  fw_r.write('#include <stdlib.h>\n')
  fw_r.write('#include <omp.h>\n')
  fw_r.write('\n')
  fw_r.write('#include "batched_blas_common.h"\n')
  fw_r.write('#include "batched_blas_fp16.h"\n')
  fw_r.write('#include "batched_blas_cost.h"\n')
  fw_r.write('#include "batched_blas_schedule.h"\n')
  fw_r.write('#include "bblas_error.h"\n')
  fw_r.write('\n')

  return

# routine declaration
def routine_part2(fw_r,row1,row2):
  routine_declaration(fw_r,row1,row2)
  fw_r.write('\n')
  fw_r.write('{\n')
  return

def routine_part3(fw_r):
  fw_r.write('  int total_batch_count;\n')
  fw_r.write('  int *which_thread;\n') 
  fw_r.write('  const int num_threads = omp_get_max_threads();\n')
  fw_r.write('  int local_no, group_no;\n')
  fw_r.write('  int my_tno;\n')
  fw_r.write('  int group_head[group_count];\n') 
  fw_r.write('  struct _cost_param cost_param;\n')
  fw_r.write('\n')
  fw_r.write('  if (group_count < 0) {\n')
  fw_r.write('    bblas_error("Illegal value of group_count");\n')
  fw_r.write('    info[0] = -1;\n')
  fw_r.write('    return;\n')
  fw_r.write('  }\n')
  fw_r.write('  int offset = 0;\n')
  fw_r.write('  int info_offset = 0;\n')
  fw_r.write('  int info_option = info[0];\n')
  fw_r.write('  int flag = 0;\n')
  fw_r.write('\n')
  fw_r.write('  total_batch_count = 0;\n')
  fw_r.write('  for (group_no = 0; group_no < group_count; group_no++) {\n')
  fw_r.write('    total_batch_count += group_size[group_no];\n') 
  fw_r.write('  }\n')
  fw_r.write('\n')
  fw_r.write('  group_head[0] = 0;\n')
  fw_r.write('  for(group_no = 1; group_no < group_count; group_no++){\n')
  fw_r.write('    group_head[group_no] = group_head[group_no - 1] + group_size[group_no - 1];\n') 
  fw_r.write('  }\n')
  fw_r.write('\n')

  return

def routine_errorcheck_part1(omp):
  fw_r.write('\n') 
  fw_r.write('      int info_local;\n') 
  fw_r.write('      if (info_option == BblasErrorsReportAll) \n') 
  fw_r.write('        info_offset = offset+1;\n') 
  fw_r.write('      else if (info_option == BblasErrorsReportGroup)\n') 
  fw_r.write('        info_offset = group_no+1;\n') 
  fw_r.write('      else\n') 
  fw_r.write('        info_offset = 0;\n') 
  fw_r.write('      info_local = info_option;\n') 
  fw_r.write('      if (group_size[group_no] < 0) {\n') 
  if omp == 'on' :
    fw_r.write('        #pragma omp master \n')
  fw_r.write('        { \n')
  fw_r.write('        bblas_error("Illegal values of group_sizes");\n') 
  fw_r.write('        info[0] = -2;\n') 
  fw_r.write('        } \n')

  return

def routine_errorcheck_part2(row2):
  if row2.count('layout') > 0 :
    fw_r.write('      if ((layout != BblasRowMajor) &&\n')
    fw_r.write('          (layout != BblasColMajor)) {\n')
    fw_r.write('          bblas_error("Illegal value of layout");\n')
    fw_r.write('          if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('            bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 1);\n')
    fw_r.write('          }\n')
    fw_r.write('        return;\n')
    fw_r.write('      }\n')

  if row2.count('transa') > 0 :
    fw_r.write('      if ((transa[group_no] != BblasNoTrans) &&\n')
    fw_r.write('	  (transa[group_no] != BblasTrans) &&\n')
    fw_r.write('	  (transa[group_no] != BblasConjTrans)) {\n')
    fw_r.write('	  bblas_error("Illegal value of transa");\n')
    fw_r.write('	  if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	    bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 2);\n')
    fw_r.write('	  }\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('transb') > 0 :
    fw_r.write('      if ((transb[group_no] != BblasNoTrans) &&\n')
    fw_r.write('	  (transb[group_no] != BblasTrans) &&\n')
    fw_r.write('	  (transb[group_no] != BblasConjTrans)) {\n')
    fw_r.write('	bblas_error("Illegal value of transb");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 3);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('m') > 0 :
    fw_r.write('      if (m[group_no] < 0) {\n')
    fw_r.write('	bblas_error("Illegal value of m");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 4);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('n') > 0 :
    fw_r.write('      if (n[group_no] < 0) {\n')
    fw_r.write('	bblas_error("Illegal value of n");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 5);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('k') > 0 :
    fw_r.write('      if (k[group_no] < 0) {\n')
    fw_r.write('	bblas_error("Illegal value of k");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 6);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('lda') > 0 :
    fw_r.write('      if (lda[group_no] < imax(1, lda[group_no])) {\n')
    fw_r.write('	bblas_error("Illegal value of lda");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 7);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('ldb') > 0 :
    fw_r.write('      if (ldb[group_no] < imax(1, ldb[group_no])) {\n')
    fw_r.write('	bblas_error("Illegal value of ldb");\n')
    fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 8);\n')
    fw_r.write('	}\n')
    fw_r.write('	return;\n')
    fw_r.write('      }\n')

  if row2.count('ldc') > 0 :
    if row2.count('m') > 0 :
      fw_r.write('      if (ldc[group_no] < imax(1, m[group_no])) {\n')
      fw_r.write('	bblas_error("Illegal value of ldc");\n')
      fw_r.write('	if (info[info_offset] != BblasErrorsReportNone) {\n')
      fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 9);\n')
      fw_r.write('	}\n')
      fw_r.write('	return;\n')
      fw_r.write('      }\n')

    return

def routine_errorcheck_part3(row2):
  if row2.count('layout') > 0 :
    fw_r.write('      if ((layout != BblasRowMajor) &&\n')
    fw_r.write('          (layout != BblasColMajor)) {\n')
    fw_r.write('          #pragma omp master \n')
    fw_r.write('          { \n')
    fw_r.write('          bblas_error("Illegal value of layout");\n')
    fw_r.write('          if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('            bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 1);\n')
    fw_r.write('          }\n')
    fw_r.write('          }\n')
    fw_r.write('        break;\n')
    fw_r.write('      }\n')

  if row2.count('transa') > 0 :
    fw_r.write('      if ((transa[group_no] != BblasNoTrans) &&\n')
    fw_r.write('	  (transa[group_no] != BblasTrans) &&\n')
    fw_r.write('	  (transa[group_no] != BblasConjTrans)) {\n')
    fw_r.write('          #pragma omp master \n')
    fw_r.write('          { \n')
    fw_r.write('	  bblas_error("Illegal value of transa");\n')
    fw_r.write('	  if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	    bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 2);\n')
    fw_r.write('	  }\n')
    fw_r.write('          }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('transb') > 0 :
    fw_r.write('      if ((transb[group_no] != BblasNoTrans) &&\n')
    fw_r.write('	  (transb[group_no] != BblasTrans) &&\n')
    fw_r.write('	  (transb[group_no] != BblasConjTrans)) {\n')
    fw_r.write('          #pragma omp master \n')
    fw_r.write('          { \n')
    fw_r.write('	  bblas_error("Illegal value of transb");\n')
    fw_r.write('	  if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	    bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 3);\n')
    fw_r.write('	  }\n')
    fw_r.write('          }\n')
    fw_r.write('	  break;\n')
    fw_r.write('        }\n')

  if row2.count('m') > 0 :
    fw_r.write('      if (m[group_no] < 0) {\n')
    fw_r.write('        #pragma omp master \n')
    fw_r.write('        { \n')
    fw_r.write('	bblas_error("Illegal value of m");\n')
    fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 4);\n')
    fw_r.write('	}\n')
    fw_r.write('        }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('n') > 0 :
    fw_r.write('      if (n[group_no] < 0) {\n')
    fw_r.write('        #pragma omp master \n')
    fw_r.write('        { \n')
    fw_r.write('	bblas_error("Illegal value of n");\n')
    fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 5);\n')
    fw_r.write('	}\n')
    fw_r.write('        }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('k') > 0 :
    fw_r.write('      if (k[group_no] < 0) {\n')
    fw_r.write('        #pragma omp master \n')
    fw_r.write('        { \n')
    fw_r.write('	bblas_error("Illegal value of k");\n')
    fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 6);\n')
    fw_r.write('	}\n')
    fw_r.write('        }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('lda') > 0 :
    fw_r.write('      if (lda[group_no] < imax(1, lda[group_no])) {\n')
    fw_r.write('        #pragma omp master \n')
    fw_r.write('        { \n')
    fw_r.write('	bblas_error("Illegal value of lda");\n')
    fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 7);\n')
    fw_r.write('	}\n')
    fw_r.write('        }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('ldb') > 0 :
    fw_r.write('      if (ldb[group_no] < imax(1, ldb[group_no])) {\n')
    fw_r.write('        #pragma omp master \n')
    fw_r.write('        { \n')
    fw_r.write('	bblas_error("Illegal value of ldb");\n')
    fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
    fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 8);\n')
    fw_r.write('	}\n')
    fw_r.write('        }\n')
    fw_r.write('	break;\n')
    fw_r.write('      }\n')

  if row2.count('ldc') > 0 :
    if row2.count('m') > 0 :
      fw_r.write('      if (ldc[group_no] < imax(1, m[group_no])) {\n')
      fw_r.write('      #pragma omp master \n')
      fw_r.write('      { \n')
      fw_r.write('	bblas_error("Illegal value of ldc");\n')
      fw_r.write('	if (info_local != BblasErrorsReportNone) {\n')
      fw_r.write('	  bblas_set_info(info[info_offset], &info[info_offset], group_size[group_no], 9);\n')
      fw_r.write('	}\n')
      fw_r.write('      }\n')
      fw_r.write('	break;\n')
      fw_r.write('      }\n')

    return


def routine_part3_1(fw_r,row1,row2):

  fw_r.write('  if(!use_batch()){\n')
  fw_r.write('    for (group_no = 0; group_no < group_count; group_no++){\n')
  routine_errorcheck_part1('off')
  fw_r.write('        return;\n') 
  fw_r.write('      }\n') 
  routine_errorcheck_part2(row2)
  fw_r.write('\n')   
  fw_r.write('      for (local_no = 0; local_no < group_size[group_no]; local_no++){\n')
  fw_r.write('        ')
  routine_call(fw_r,row1,row2)
  fw_r.write('        if (info_local == BblasErrorsReportAll)\n')
  fw_r.write('          info[local_no+info_offset] = 0;\n')
  fw_r.write('      }\n')
  fw_r.write('      if (info_local != BblasErrorsReportAll)\n')
  fw_r.write('        info_local = 0;\n')
  fw_r.write('      if (info_local != 0 && flag == 0) {\n')
  fw_r.write('        info[0] = info[info_offset];\n')
  fw_r.write('        flag = 1;\n')
  fw_r.write('      }\n')
  fw_r.write('      offset += group_size[group_no];\n')
  fw_r.write('    }\n')
  fw_r.write('    return;\n')
  fw_r.write('  }\n')
  fw_r.write('\n')

  return

def routine_part4(fw_r,row3):
  
  c = 0
  for r in row3:
    if c == 0:
      pass
    else:
      fw_r.write('  cost_param.ptr[')
      fw_r.write('%d' % (c-1))
      fw_r.write('] = (void *)')
      fw_r.write(r)
      fw_r.write(';')
      fw_r.write('\n')
    c = c + 1
    
  fw_r.write('\n')
  fw_r.write('\n')

  return

def routine_part5(fw_r,row3):
  fw_r.write('  which_thread = (int *) malloc(sizeof(int) * total_batch_count);\n')
  fw_r.write('  schedule_batch(group_count, group_size,')
  fw_r.write(row3[0])
  fw_r.write(',&cost_param, which_thread);\n')
  fw_r.write('\n')
  fw_r.write('#pragma omp parallel default(shared) num_threads(num_threads)	\\\n')
  fw_r.write('  private(group_no, local_no, my_tno) firstprivate(flag, offset, info_offset) \n')
  fw_r.write('  {\n')
  fw_r.write('  my_tno = omp_get_thread_num();\n')
  fw_r.write('  for (group_no = 0; group_no < group_count; group_no++){\n')
  routine_errorcheck_part1('on')
  fw_r.write('        break;\n') 
  fw_r.write('      }\n') 
  routine_errorcheck_part3(row2)
  fw_r.write('\n')   
  fw_r.write('    for (local_no = 0; local_no < group_size[group_no]; local_no++){\n')
  fw_r.write('      if(which_thread[group_head[group_no]+local_no] == my_tno){\n')
  
  return

# subroutine call
def routine_part6(fw_r,row1,row2):
  fw_r.write('        ')
  
  routine_call(fw_r,row1,row2)
  
  return

def routine_part7(fw_r):
  fw_r.write('        if (info_local == BblasErrorsReportAll)\n')
  fw_r.write('          info[local_no+info_offset] = 0;\n')
  fw_r.write('      }\n')
  fw_r.write('    }\n')
  fw_r.write('      if (info_local != BblasErrorsReportAll)\n')
  fw_r.write('        info_local = 0;\n')
  fw_r.write('      if (info_local != 0 && flag == 0) {\n')
  fw_r.write('        #pragma omp critical \n')
  fw_r.write('        { \n')
  fw_r.write('        info[0] = info[info_offset];\n')
  fw_r.write('        } \n')
  fw_r.write('        flag = 1;\n')
  fw_r.write('      }\n')
  fw_r.write('  }\n')
  fw_r.write('      offset += group_size[group_no];\n')
  fw_r.write('  }\n')
  fw_r.write('\n')
  fw_r.write('  free(which_thread);\n')
  fw_r.write('\n')
  fw_r.write('  return;\n')
  fw_r.write('}\n')
  fw_r.write('\n')

  return

def header_part1(fw_h):
  fw_h.write('#ifndef _BATCHED_BLAS_APRIORICOST_H_\n')
  fw_h.write('#define _BATCHED_BLAS_APRIORICOST_H_\n')
  fw_h.write('\n')
  fw_h.write('#include "batched_blas_common.h"\n')
  fw_h.write('\n')

  return

def header_part2(fw_h,row1,row2):
  routine_declaration(fw_h,row1,row2)
  fw_h.write(';')
  fw_h.write('\n')
  return

def header_part3(fw_h):
  fw_h.write('\n')
  fw_h.write('#endif\n')

  return

def make_part1(fw_mk):
  fw_mk.write('ifeq (x$(CC),xcc)\n')
  fw_mk.write('ifdef MKLROOT\n')
  fw_mk.write('CC=icc\n')
  fw_mk.write('else\n')
  fw_mk.write('ifdef FJSVXTCLANGA\n')
  fw_mk.write('CC=fccpx\n')
  fw_mk.write('else\n')
  fw_mk.write('CC=gcc\n')
  fw_mk.write('endif\n')
  fw_mk.write('endif\n')
  fw_mk.write('endif\n')
  fw_mk.write('\n')
  fw_mk.write('ifeq (x$(BASEBLAS),xMKL)\n')
  fw_mk.write('CCFLAG=-O3 -xHost -qopenmp -I./\n')
  fw_mk.write('BLAS=-lmkl_rt -lm\n')
  fw_mk.write('else\n')
  fw_mk.write('ifeq (x$(BASEBLAS),xSSL2)\n')
  fw_mk.write('CCFLAG=-std=gnu11 -Kfast,ocl,openmp -Nclang -D_CBLAS_ -I./\n')
  fw_mk.write('BLAS=-SSL2BLAMP -lm\n')
  fw_mk.write('else\n')
  fw_mk.write('ifeq (x$(BASEBLAS),xARMPL)\n')
  fw_mk.write('ARMPL_DIR=/opt/arm/armpl-21.0.0_AArch64-SVE_RHEL-7_gcc_aarch64-linux\n')
  fw_mk.write('CCFLAG=-O3 -fopenmp -Wall -D_CBLAS_ -I./ -I${ARMPL_DIR}/include -D__ARM_ACLE -D_ARMPL_\n')
  fw_mk.write('BLAS=${ARMPL_DIR}/lib/libarmpl_lp64.a -lgfortran -lm\n')
  fw_mk.write('else\n')
  fw_mk.write('CCFLAG=-O3 -fopenmp -Wall -D_CBLAS_ -I./ -I/usr/include/openblas\n')
  fw_mk.write('BLAS=-lopenblas -lm\n')
  fw_mk.write('endif\n')
  fw_mk.write('endif\n')
  fw_mk.write('endif\n')
  fw_mk.write('\n')
  fw_mk.write('CCFLAG+=-fPIC\n')
  fw_mk.write('\n')
  #fw_mk.write('CC=icc\n')
  #fw_mk.write('CCFLAG=-O3 -xHost -qopenmp -I./ \n')
  #fw_mk.write('#CC=gcc\n')
  #fw_mk.write('#CCFLAG=-O3 -fopenmp -Wall -D_CBLAS_ -I./ \n')
  #fw_mk.write('#CC=fccpx\n')
  #fw_mk.write('#CCFLAG=-Kfast,ocl,openmp -Nclang -D_CBLAS_ -I./ \n')
  fw_mk.write('\n')
  fw_mk.write('AR=ar\n')
  #fw_mk.write('BLAS=-mkl=parallel\n')
  fw_mk.write('\n')
  #fw_mk.write('LIB_BATCHED_BLAS = libbatched_blas.a\n')
  fw_mk.write('LIB_BATCHED_BLAS = libbblas_aprioricost.a\n')
  fw_mk.write('LIB_BATCHED_BLAS_SO = libbblas_aprioricost.so\n')
  fw_mk.write('\n')
  fw_mk.write('OBJS_COMMON = \\\n')
  fw_mk.write('batched_blas_consts.o \\\n')
  fw_mk.write('batched_blas_common.o \\\n')
  fw_mk.write('batched_blas_schedule.o \\\n')
  fw_mk.write('batched_blas_cost.o\n')
  fw_mk.write('\n')
  fw_mk.write('OBJS_SRC = ')

  return

def make_part2(fw_mk,routine_name):
  fw_mk.write(' \\\n')
  fw_mk.write(routine_name + '.o')
  return

def make_part3(fw_mk):
  fw_mk.write('\n')
  fw_mk.write('\n')
  fw_mk.write('OBJS = $(OBJS_COMMON) $(OBJS_SRC)\n')
  fw_mk.write('\n')
  fw_mk.write('all : $(LIB_BATCHED_BLAS) $(LIB_BATCHED_BLAS_SO)\n')
  fw_mk.write('$(LIB_BATCHED_BLAS) : $(OBJS)\n')
  fw_mk.write('\t$(AR) r $(LIB_BATCHED_BLAS) $(OBJS)\n')
  fw_mk.write('$(LIB_BATCHED_BLAS_SO) : $(OBJS)\n')
  fw_mk.write('\t$(CC) -fPIC -shared -o $(LIB_BATCHED_BLAS_SO) $(OBJS)\n')
  fw_mk.write('\n')
  fw_mk.write('.c.o :\n')
  fw_mk.write('\t$(CC) $(CCFLAG) -c $< $(DEF)\n')
  fw_mk.write('\n')
  fw_mk.write('clean :\n')
  fw_mk.write('\trm -f $(LIB_BATCHED_BLAS) $(LIB_BATCHED_BLAS_SO) *.o\n')

  return

def copy_files():
  shutil.copy('include/batched_blas_common.c', src_path)
  shutil.copy('include/batched_blas_common.h', src_path)
  shutil.copy('include/batched_blas_cost.c', src_path)
  shutil.copy('include/batched_blas_cost.h', src_path)
  shutil.copy('include/batched_blas_schedule.c', src_path)
  shutil.copy('include/batched_blas_schedule.h', src_path)
  shutil.copy('include/bblas_aprioricost.h', src_path)
  shutil.copy('include/bblas_types.h', src_path)
  shutil.copy('include/bblas_error.h', src_path)
  shutil.copy('include/batched_blas_consts.c', src_path)
  shutil.copy('include/batched_blas_fp16.h', src_path)
  return

def routne_cblas_param(iofile):
  ofile = iofile
  ifile = iofile + '_'
  shutil.copyfile(iofile,ifile)
  fr = open(ifile,'r')
  fw = open(ofile,'w')

  for row in fr:
    if 'cblas_' in row:
      row_rep = row
      if 'diag[group_no]' in row_rep:
        row_rep = row_rep.replace('diag[group_no]','diag_cblas( diag[group_no] ) ')
      if 'side[group_no]' in row_rep:                
        row_rep = row_rep.replace('side[group_no]','side_cblas( side[group_no] ) ')
      if 'trans[group_no]' in row_rep:                
        row_rep = row_rep.replace('trans[group_no]','transpose_cblas( trans[group_no] ) ')
      if 'transa[group_no]' in row_rep:                
        row_rep = row_rep.replace('transa[group_no]','transpose_cblas( transa[group_no] ) ')
      if 'transb[group_no]' in row_rep:                
        row_rep = row_rep.replace('transb[group_no]','transpose_cblas( transb[group_no] ) ')
      if 'uplo[group_no]' in row_rep:                
        row_rep = row_rep.replace('uplo[group_no]','uplo_cblas( uplo[group_no] ) ')
      if 'layout' in row_rep:                
        row_rep = row_rep.replace('layout','layout_cblas( layout ) ')
      fw.write(row_rep)
    else:
      fw.write(row)

  fr.close()
  fw.close()
  os.remove(ifile)
  return

args = sys.argv
argc = len(args)

if (argc < 2):
  print('Usage: python *.py <description_file.csv>')
  quit()

fr = open(args[1],'r')

reader = csv.reader(fr)

# create header file
fw_h = open(src_path + '/' + 'batched_blas_aprioricost.h','w')
# create Make file
fw_mk = open(src_path + '/' + 'Makefile','w')

header_part1(fw_h)
make_part1(fw_mk)

mode = 0
for row in reader:
  if mode == 0:
    row1 = row
  elif mode == 1:
    row2 = row
  else:
    row3 = row
    routine_type = row1[1]
    routine_name = row1[2]
    # make c program
    ofile = src_path + '/' + row1[0] + '.c'
    fw_r = open(ofile,'w')

    if 'my_sgemm' in row1[0]:
      if '_batchf' in row1[0] :
        routine_part1(fw_r)
        fw_r.write('#include "batched_blas_aprioricost.h"\n')
        routine_part2(fw_r,row1,row2)
        routine_fixed(fw_r,row1,row2)

      else :
        routine_part1(fw_r)
        routine_my_sgemm(fw_r)
        routine_part2(fw_r,row1,row2)
        routine_part3(fw_r)
        routine_part3_1(fw_r,row1,row2)
        routine_part4(fw_r,row3)
        routine_part5(fw_r,row3)
        routine_part6(fw_r,row1,row2)
        routine_part7(fw_r)
    else :
      if '_batchf' in row1[0] :
        routine_part1(fw_r)
        fw_r.write('#include "batched_blas_aprioricost.h"\n')
        routine_part2(fw_r,row1,row2)
        routine_fixed(fw_r,row1,row2)

      else :
        routine_part1(fw_r)
        routine_part2(fw_r,row1,row2)
        routine_part3(fw_r)
        routine_part3_1(fw_r,row1,row2)
        routine_part4(fw_r,row3)
        routine_part5(fw_r,row3)
        routine_part6(fw_r,row1,row2)
        routine_part7(fw_r)

    fw_r.close()
    routne_cblas_param(ofile)
    # make header file
    header_part2(fw_h,row1,row2)
    # make make file
    make_part2(fw_mk,row1[0])
  
  mode = (mode + 1) % 3
  
make_part3(fw_mk)
fw_mk.close()

header_part3(fw_h)
fw_h.close()

fr.close()

copy_files()
