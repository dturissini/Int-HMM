from __future__ import division
from datetime import datetime
import random
import math
import sys
import MySQLdb
import os
from warnings import filterwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

cnx = MySQLdb.connect(user='', passwd='', db='', host='')



#model variable names were chosen to match those in chapter 3 of "Biological sequence analysis" by Durbin et al.

def forward():
  for i in range(1, len(x)):
    f_tilde_sums = [{} for z in range(0, len(x))]
    s[i] = 0
    for l in states:
      f_tilde_sums[i][l] = 0        
      if i == 1:    
        for k in a[i].keys():
          f_tilde_sums[i][l] += f_tilde[i - 1][k] * a[i][k][l]
      else:
        for k in states:
          f_tilde_sums[i][l] += f_tilde[i - 1][k] * a[i][k][l]        
      s[i] += e[i][l] * f_tilde_sums[i][l]  
                  
    for l in states:    
      f_tilde[i][l] = e[i][l] * f_tilde_sums[i][l] / s[i]  


def backward():
  for i in range(len(x) - 2, 0, -1):
    for k in states:
      b_tilde[i][k] = 0
      for l in states:
        b_tilde[i][k] += a[i + 1][k][l] * b_tilde[i + 1][l] * e[i + 1][l] / s[i]
      

def nCr(n,k):
    return round(math.exp(math.log(math.factorial(n)) - math.log(math.factorial(k)) - math.log(math.factorial(n-k))))


def binomial(n, k, p):
  return nCr(n,k) * p ** (k) * (1 - p) ** (n - k)


def allele_prob(n, k, p, pa, pb):
  prob = 0
  for i in range(0, (n + 1)):
    pab_sum = 0
    for j in range(max(0, k - i), min(k, n - i) + 1):
      pab_sum += nCr(n - i, j) * pa ** j * (1 - pa) ** (n - i - j) * nCr(i, j - k + i) * pb ** (j - k + i) * (1 - pb) ** (k - j)
    prob += nCr(n, i) * p ** i * (1 - p) ** (n - i) * pab_sum
  return prob



#start
(ind, prefix, arm, pab, base_a_prob, err_multiplier, freq_diff_cutoff) = sys.argv[1:]

print "\t".join([datetime.now().strftime("%y-%m-%d %H:%M:%S"), 'start'])



pab = float(pab)
base_a_prob = float(base_a_prob)
hmm_table = prefix + '_' + ind.replace(".", "_") + '_hmm_dsm_pab_' + '0' * (2 - len(str(int(pab * 100)))) + str(int(pab * 100)) + '_' + str(int(-1 * math.log(base_a_prob, 10) + .5)) + '_' + err_multiplier + '_all_' + str(int(float(freq_diff_cutoff) * 100))

if pab < .01:
  hmm_table = prefix + '_' + ind.replace(".", "_") + '_hmm_dsm_pab_' + '0' * (3 - len(str(int(pab * 1000)))) + str(int(pab * 1000)) + '_' + str(int(-1 * math.log(base_a_prob, 10) + .5)) + '_' + err_multiplier + '_all_' + str(int(float(freq_diff_cutoff) * 100))

err_multiplier = float(err_multiplier)

sex = 'female'

max_cov = 25

                            
genotype_table = 'sites_all_' + prefix + '_' + str(int(float(freq_diff_cutoff) * 100)) + '_' + arm

genotype_cursor = cnx.cursor()
genotype_cursor.execute("""select pos, ref_allele_cov + alt_allele_cov, ref_allele_cov, alt_allele_cov, recip_ref_freq, donor_ref_freq
                           from %s
                           where ind = '%s'
                           order by pos""" % (genotype_table, ind))     

print(arm)
positions = [0] 
x = [None]
for (pos, cov, a_cov, b_cov, par_1_freq, par_2_freq) in genotype_cursor:
  if cov > max_cov:
    a_cov = int(round(max_cov * (a_cov / cov)))
    b_cov = max_cov - a_cov
    cov = max_cov
  x.append([cov, a_cov, b_cov, float(par_1_freq), float(par_2_freq)])
  positions.append(pos)


pa = pab
pb = pab  
e = [{} for z in range(0, len(x))]
states = ['homo_1', 'het', 'homo_2', 'err_homo_1', 'err_het', 'err_homo_2']
for i in range(1,len(x)):
  for state_k in states: 
    if arm == 'X' and sex == 'male':
      if state_k == 'homo_1' or state_k == 'err_homo_1':
        e[i][state_k] = allele_prob(x[i][0], x[i][1], x[i][3], pa, pb)
      elif state_k == 'homo_2' or state_k == 'err_homo_2':
        e[i][state_k] = allele_prob(x[i][0], x[i][1], x[i][4], pa, pb)
      else:
        e[i][state_k] = 0        
    else:
      if state_k == 'homo_1' or state_k == 'err_homo_1':
        e[i][state_k] = allele_prob(x[i][0], x[i][1], x[i][3], pa, pb)
      elif state_k == 'homo_2' or state_k == 'err_homo_2':
        e[i][state_k] = allele_prob(x[i][0], x[i][1], x[i][4], pa, pb)
      else:
        e[i][state_k] = allele_prob(x[i][0], x[i][1], (x[i][3] + x[i][4]) / 2, pa, pb)


a = [{} for z in range(0, len(x))]
a[1] = {0: {}}
for state_1 in states:
  a[1][0][state_1] = 1 / len(states)
  for i in range(1,len(positions)):
    poisson_lambda = base_a_prob * float(abs(positions[i - 1] - positions[i]))
    a_prob = poisson_lambda * math.exp(-poisson_lambda)
    
    poisson_lambda_err = err_multiplier * base_a_prob * float(abs(positions[i - 1] - positions[i]))
    a_prob_err = poisson_lambda_err * math.exp(-poisson_lambda_err)
    
    a[i][state_1] = {}
    if state_1 == 'homo_1':
      a[i][state_1]['homo_1'] = 1 - 4 * a_prob
      a[i][state_1]['het'] = a_prob
      a[i][state_1]['homo_2'] = a_prob
      a[i][state_1]['err_homo_1'] = 0
      a[i][state_1]['err_het'] = a_prob
      a[i][state_1]['err_homo_2'] = a_prob
    elif state_1 == 'het':
      a[i][state_1]['homo_1'] = a_prob
      a[i][state_1]['het'] = 1 - 13 * a_prob
      a[i][state_1]['homo_2'] = a_prob
      a[i][state_1]['err_homo_1'] = 10 * a_prob
      a[i][state_1]['err_het'] = 0
      a[i][state_1]['err_homo_2'] = a_prob
    elif state_1 == 'homo_2':
      a[i][state_1]['homo_1'] = a_prob
      a[i][state_1]['het'] = a_prob
      a[i][state_1]['homo_2'] = 1 - 13 * a_prob
      a[i][state_1]['err_homo_1'] = 10 * a_prob
      a[i][state_1]['err_het'] = a_prob
      a[i][state_1]['err_homo_2'] = 0
    elif state_1 == 'err_homo_1':
      a[i][state_1]['homo_1'] = 0
      a[i][state_1]['het'] = (1 - 2 * a_prob_err) / 4
      a[i][state_1]['homo_2'] = (1 - 2 * a_prob_err) / 4
      a[i][state_1]['err_homo_1'] = 2 * a_prob_err
      a[i][state_1]['err_het'] = (1 - 2 * a_prob_err) / 4
      a[i][state_1]['err_homo_2'] = (1 - 2 * a_prob_err) / 4
    elif state_1 == 'err_het':
      a[i][state_1]['homo_1'] = (1 - a_prob_err) / 4
      a[i][state_1]['het'] = 0
      a[i][state_1]['homo_2'] = (1 - a_prob_err) / 4
      a[i][state_1]['err_homo_1'] = 9 / 10 * (1 - a_prob_err) / 2
      a[i][state_1]['err_het'] = a_prob_err
      a[i][state_1]['err_homo_2'] = 1 / 10 * (1 - a_prob_err) / 2
    elif state_1 == 'err_homo_2':
      a[i][state_1]['homo_1'] = (1 - a_prob_err) / 4
      a[i][state_1]['het'] = (1 - a_prob_err) / 4
      a[i][state_1]['homo_2'] = 0
      a[i][state_1]['err_homo_1'] = 9 / 10 * (1 - a_prob_err) / 2
      a[i][state_1]['err_het'] = 1 / 10 * (1 - a_prob_err) / 2
      a[i][state_1]['err_homo_2'] = a_prob_err


s = [None] * (len(x))

f = [{} for z in range(0, len(x))]
f_tilde = [{} for z in range(0, len(x))]
b = [{} for z in range(0, len(x))]
b_tilde = [{} for z in range(0, len(x))]

f[0] = {0: 1}
f_tilde[0] = {0: 1}
for state in states:
  f[0][state] = 0
  f_tilde[0][state] = 0 
  b[len(x) - 1][state] = 1


forward()


for k in states:
  b_tilde[len(x) - 1][k] = b[len(x) - 1][k] / s[len(x) - 1]
    
    
backward()


f_tilde_L_sum = 0
for k in states:
  f_tilde_L_sum += f_tilde[len(x) - 1][k]


log_s_prod_f = 0
log_s_prod_b = 0
for i in range (1, len(x)):
  log_s_prod_b += math.log(s[i])

cnx = MySQLdb.connect(user='matute_user', passwd='drosophila', db='matute_db', host='bioapps.its.unc.edu')
insert_cursor = cnx.cursor()
log_s_prod_f_L = log_s_prod_b
Pkx = [{} for z in range(0, len(x))]
for i in range(1, len(x)):
  log_s_prod_f += math.log(s[i])
  if i > 1:
    log_s_prod_b -= math.log(s[i - 1])
    
  for k in states:
    if f_tilde[i][k] > 0 and b_tilde[i][k] > 0:
      Pkx[i][k] = math.exp(math.log(f_tilde[i][k]) + math.log(b_tilde[i][k]) + log_s_prod_f + log_s_prod_b - log_s_prod_f_L - math.log(f_tilde_L_sum))  
    else:
      Pkx[i][k] = 0
  
  
  probs = (Pkx[i]['homo_1'], Pkx[i]['het'], Pkx[i]['homo_2'], Pkx[i]['err_homo_1'], Pkx[i]['err_het'], Pkx[i]['err_homo_2'])
  prob_names = ('homo_1', 'het', 'homo_2', 'err_homo_1', 'err_het', 'err_homo_2')
  max_prob = max(probs)
  max_prob_index = probs.index(max_prob)
  
  insert_cursor.execute("""insert into %s
                           values
                           (null, '%s', '%s', %s, %s, %s, %s, %s, %s, %s, %s, %s, '%s')""" % (hmm_table, ind, arm, positions[i], Pkx[i]['homo_1'], Pkx[i]['het'], Pkx[i]['homo_2'], Pkx[i]['err_homo_1'], Pkx[i]['err_het'], Pkx[i]['err_homo_2'], x[i][1], x[i][2], prob_names[max_prob_index]))
  cnx.commit()                         
print("done")
