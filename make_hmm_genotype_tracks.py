from __future__ import division
import MySQLdb
import sys
from datetime import datetime
from warnings import filterwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

cnx = MySQLdb.connect(user='', passwd='', db='', host='')


def filter_tracks(filter_name, tmp_table, new_tracks_table, arm, ind, filtered_ict_ids, hmm_table): 
  filter_log_table = hmm_table + '_fillog'
  
  track_cursor = cnx.cursor()
  track_cursor.execute("""select ict_id 
                          from %s
                          where arm = '%s'
                          and ind = '%s'
                          order by ict_id""" % (tmp_table, arm, ind))


  for (ict_id,) in track_cursor:
    filtered = 'N'
    for ict_id_start in filtered_ict_ids:
      if ict_id >= ict_id_start and ict_id <= filtered_ict_ids[ict_id_start]['end']:
        filtered = 'Y'
        if ict_id == ict_id_start:
          combine_cursor = cnx.cursor()
          combine_cursor.execute("""insert into %s 
                                    select null, ind, '%s',
                                           arm, min(start), max(end),
                                           sum(total_snps), sum(introgress_snps), max(end) - min(start) + 1, sum(cum_cov), 0
                                    from %s
                                    where arm = '%s'
                                    and ind = '%s'
                                    and ict_id between %s and %s
                                    group by ind, arm, '%s'""" % (new_tracks_table, filtered_ict_ids[ict_id_start]['genotype'], 
                                                                  tmp_table, arm, ind, 
                                                                  ict_id_start, filtered_ict_ids[ict_id_start]['end'],
                                                                  filtered_ict_ids[ict_id_start]['genotype']))
                                                                      
          cnx.commit()
          
          num_tracks = filtered_ict_ids[ict_id_start]['end'] - ict_id_start + 1
          log_cursor = cnx.cursor()
          log_cursor.execute("""insert into %s
                                values
                                (null, '%s', '%s', %s, '%s')""" % (filter_log_table, ind, filter_name, num_tracks, filtered_ict_ids[ict_id_start]['genotype']))
          cnx.commit()



    if filtered == 'N':
      insert_track_cursor = cnx.cursor()
      insert_track_cursor.execute("""insert into %s 
                                     select null, ind, genotype, arm, start, end, total_snps, introgress_snps, track_len, cum_cov, per_repeat
                                     from %s
                                     where ict_id = %s""" % (new_tracks_table, tmp_table, ict_id))
      cnx.commit() 
                                    
    if ict_id in filtered_ict_ids.values():
      for ict_id_start in filtered_ict_ids:
        if ict_id == filtered_ict_ids[ict_id_start]['end']:
          filtered_ict_ids.pop(ict_id_start)
          break
  

                                  

def merge_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmptrk'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         per_repeat decimal(5,2),
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 


  merge_cursor = cnx.cursor()
  #note that t.genotype is included in the group by solely to get the mysql optimizer to range check on indexes for t2 and t3
  #also the code as is may not always produce the desired genotype when adjacent pairs are merged
  #100(homo2) + 5(homo2) + 6(het) -> het
  #this will need to be addresed
  merge_cursor.execute("""select straight_join t.ict_id, t3.ict_id, if(if(t.genotype = 'homo_2', t.total_snps, 0) + if(t3.genotype = 'homo_2', t3.total_snps, 0) >= if(t.genotype = 'het', t.total_snps, 0) + if(t3.genotype = 'het', t3.total_snps, 0), 'homo_2', 'het'), t.genotype
                          from %s t, %s t3, %s t2
                          where t3.ict_id between t.ict_id + 2 and t.ict_id + 50
                          and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                          and t.genotype in ('het', 'homo_2')
                          and t3.genotype in ('het', 'homo_2')
                          group by t.ict_id, t3.ict_id, t.genotype, t.total_snps, t3.genotype, t3.total_snps
                          having sum(if(t2.genotype like 'err%%', 1, 0)) = count(*)
                          and sum(if(t2.genotype = 'err_homo_1', t2.total_snps, 0)) <= 1 + sum(if(t2.genotype in ('err_het', 'err_homo_2'), t2.total_snps, 0))
                          order by t.ict_id
                          """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype, ignore_genotype) in merge_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  filter_tracks('merge', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  




  
  
  
def remove_err_filter(tracks_table, ind, arm, hmm_table, search_size):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmptrk'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         per_repeat decimal(5,2),
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 
   
  #the relative order of the queries in the union is important and needs to be maintained here
  remove_err_cursor = cnx.cursor()
  remove_err_cursor.execute("""select straight_join t.ict_id ict_id_start, if(t3.genotype = 'homo_1', t3.ict_id, t3.ict_id - 1) ict_id_end, 'homo_1'
                               from %s t, %s t3, %s t2
                               where t3.ict_id between t.ict_id + 1 and t.ict_id + %s
                               and t2.ict_id between t.ict_id and t3.ict_id - 1
                               and t3.genotype not like 'err%%'
                               and t.genotype like 'err%%'
                               and t.ict_id in (select min(ict_id) from %s)
                               group by t.ict_id, t3.ict_id, t3.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               union all
                               select straight_join if(t.genotype = 'homo_1', t.ict_id, t.ict_id + 1), t3.ict_id, 'homo_1'
                               from %s t3, %s t, %s t2
                               where t3.ict_id in (select max(ict_id) from %s)
                               and t.ict_id between t3.ict_id - %s - 1 and t3.ict_id -1
                               and t2.ict_id between t.ict_id + 1 and t3.ict_id
                               and t.genotype not like 'err%%'
                               and t3.genotype like 'err%%'
                               group by t.ict_id, t3.ict_id, t.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               union all
                               select straight_join t.ict_id, t3.ict_id, 'homo_1'
                               from %s t, %s t3, %s t2
                               where t3.ict_id between t.ict_id + 2 and t.ict_id + %s
                               and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                               and t.genotype = 'homo_1'
                               and t.genotype = t3.genotype
                               group by t.ict_id, t3.ict_id, t.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               union all
                               select straight_join t.ict_id + 1, t3.ict_id, 'homo_1'
                               from %s t, %s t3, %s t2
                               where t3.ict_id between t.ict_id + 2 and t.ict_id + %s
                               and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                               and t.genotype in ('het', 'homo_2')
                               and t3.genotype = 'homo_1'
                               group by t.ict_id, t3.ict_id, t3.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               union all
                               select straight_join t.ict_id, t3.ict_id - 1, 'homo_1'
                               from %s t, %s t3, %s t2
                               where t3.ict_id between t.ict_id + 2 and t.ict_id + %s
                               and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                               and t.genotype = 'homo_1'
                               and t3.genotype in ('het', 'homo_2')
                               group by t.ict_id, t3.ict_id, t.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               union all
                               select straight_join t.ict_id + 1, t3.ict_id - 1, 'homo_1'
                               from %s t, %s t3, %s t2
                               where t3.ict_id between t.ict_id + 2 and t.ict_id + %s
                               and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                               and t.genotype in ('het', 'homo_2')
                               and t3.genotype in ('het', 'homo_2')
                               group by t.ict_id, t3.ict_id, t.genotype
                               having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                               order by ict_id_start, ict_id_end
                               """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table, search_size, tracks_tmp_filtered_table, 
                                      tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table, 
                                      tracks_tmp_filtered_table, search_size, tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table, 
                                      search_size, tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table,
                                      search_size, tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table,
                                      search_size, tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table, search_size))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype) in remove_err_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end

  filter_tracks('remove_err', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  


def merge_het_homo_2_tmp_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmptrk'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         per_repeat decimal(5,2),
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 
   

  merge_het_homo_1_cursor = cnx.cursor()
  #note that t.genotype is included in the group by solely to get the mysql optimizer to range check on indexes for t2 and t3
  merge_het_homo_1_cursor.execute("""select straight_join t.ict_id + 1, t3.ict_id - 1, if(sum(if(t2.genotype = 'homo_2', t2.introgress_snps, 0)) > sum(if(t2.genotype = 'het', t2.introgress_snps, 0)), 'homo_2', 'het'), t.genotype
                                     from %s t, %s t3, %s t2
                                     where t3.ict_id between t.ict_id + 3 and t.ict_id + 50
                                     and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                                     and t.genotype in ('homo_1', 'err')
                                     and t3.genotype in ('homo_1', 'err')
                                     group by t.ict_id + 1, t3.ict_id - 1, t.genotype
                                     having sum(if(t2.genotype = 'homo_2', 1, 0)) > 0
                                     and sum(if(t2.genotype = 'het', 1, 0)) > 0
                                     and sum(if(t2.genotype = 'homo_1', 1, 0)) = 0
                                     and sum(if(t2.genotype like 'err%%', 1, 0)) = 0
                                     order by t.ict_id
                                     """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype, ignore_genotype) in merge_het_homo_1_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  filter_tracks('merge_het_homo_2_tmp', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  




def merge_het_homo_2_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmptrk'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         per_repeat decimal(5,2),
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 
   

  merge_het_homo_2_cursor = cnx.cursor()
  #note that t.genotype is included in the group by solely to get the mysql optimizer to range check on indexes for t2 and t3
  merge_het_homo_2_cursor.execute("""select straight_join t.ict_id + 1, t3.ict_id - 1, if(sum(if(t2.genotype = 'homo_2', t2.introgress_snps, 0)) > sum(if(t2.genotype = 'het', t2.introgress_snps, 0)), 'homo_2', 'het'), t.genotype
                                     from %s t, %s t3, %s t2
                                     where t3.ict_id between t.ict_id + 3 and t.ict_id + 50
                                     and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                                     and t.genotype = t3.genotype
                                     and t.genotype = 'homo_1'
                                     group by t.ict_id + 1, t3.ict_id - 1, t.genotype
                                     having sum(if(t2.genotype = 'homo_2', 1, 0)) > 0
                                     and sum(if(t2.genotype = 'het', 1, 0)) > 0
                                     and sum(if(t2.genotype = 'homo_1', 1, 0)) = 0
                                     and sum(if(t2.genotype like 'err%%', 1, 0)) = 0
                                     order by t.ict_id
                                     """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype, ignore_genotype) in merge_het_homo_2_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  filter_tracks('merge_het_homo_2', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  



def merge_small_homo_1_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmptrk'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         per_repeat decimal(5,2),
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 


  merge_cursor = cnx.cursor()
  merge_cursor.execute("""select straight_join t.ict_id, t3.ict_id, if(if(t.genotype = 'homo_2', t.total_snps, 0) + if(t3.genotype = 'homo_2', t3.total_snps, 0) >= if(t.genotype = 'het', t.total_snps, 0) + if(t3.genotype = 'het', t3.total_snps, 0), 'homo_2', 'het')
                          from %s t, %s t3, %s t2
                          where t3.ict_id = t.ict_id + 2
                          and t2.ict_id = t.ict_id + 1
                          and t.genotype in ('het', 'homo_2')
                          and t2.genotype = 'homo_1'
                          and t3.genotype in ('het', 'homo_2')
                          and t.total_snps > t2.total_snps
                          and t3.total_snps > t2.total_snps
                          and t2.total_snps < 5
                          and (t.total_snps >= 10
                            or t3.total_snps >= 10)
                          order by t.ict_id
                          """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype) in merge_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  filter_tracks('merge_small_homo_1_filter', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  





def err_filter_tracks(filter_name, tmp_table, new_tracks_table, arm, ind, filtered_ict_ids, hmm_table):   
  track_cursor = cnx.cursor()
  track_cursor.execute("""select ict_id 
                          from %s
                          where arm = '%s'
                          and ind = '%s'
                          order by ict_id""" % (tmp_table, arm, ind))


  for (ict_id,) in track_cursor:
    filtered = 'N'
    for ict_id_start in filtered_ict_ids:
      if ict_id >= ict_id_start and ict_id <= filtered_ict_ids[ict_id_start]['end']:
        filtered = 'Y'
        if ict_id == ict_id_start:
          combine_cursor = cnx.cursor()
          combine_cursor.execute("""insert into %s 
                                    select null, ind, '%s',
                                           arm, min(start), max(end),
                                           sum(total_snps), sum(introgress_snps), max(end) - min(start) + 1, sum(cum_cov), sum(total_tracks)
                                    from %s
                                    where arm = '%s'
                                    and ind = '%s'
                                    and ict_id between %s and %s
                                    group by ind, arm, '%s'""" % (new_tracks_table, filtered_ict_ids[ict_id_start]['genotype'], 
                                                                  tmp_table, arm, ind, 
                                                                  ict_id_start, filtered_ict_ids[ict_id_start]['end'],
                                                                  filtered_ict_ids[ict_id_start]['genotype']))
                                                                      
          cnx.commit()
          
          num_tracks = filtered_ict_ids[ict_id_start]['end'] - ict_id_start + 1



    if filtered == 'N':
      insert_track_cursor = cnx.cursor()
      insert_track_cursor.execute("""insert into %s 
                                     select null, ind, genotype, arm, start, end, total_snps, introgress_snps, track_len, cum_cov, total_tracks
                                     from %s
                                     where ict_id = %s""" % (new_tracks_table, tmp_table, ict_id))
      cnx.commit() 
                                    
    if ict_id in filtered_ict_ids.values():
      for ict_id_start in filtered_ict_ids:
        if ict_id == filtered_ict_ids[ict_id_start]['end']:
          filtered_ict_ids.pop(ict_id_start)
          break
  



def singleton_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmpfil'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         total_tracks int,
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 
   
  singleton_cursor = cnx.cursor()
  singleton_cursor.execute("""select straight_join t.ict_id, t3.ict_id, t3.genotype
                              from %s t, %s t3, %s t2
                              where t3.ict_id between t.ict_id + 2 and t.ict_id + 52
                              and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                              and t3.genotype = 'homo_1'
                              and t.genotype = 'homo_1'
                              group by t.ict_id, t3.ict_id, t3.genotype
                              having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                              """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype) in singleton_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  err_filter_tracks('singleton', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  


                            
def merge_err_filter(tracks_table, ind, arm, hmm_table):
  tracks_tmp_filtered_table = tracks_table[:-6] + 'tmpfil'
  
  tmp_cur = cnx.cursor()
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  tmp_cur.execute("""create table %s
                        (ict_id int primary key,
                         ind varchar(30),
                         genotype varchar(20),
                         arm varchar(5),
                         start int,
                         end int,
                         total_snps int,
                         introgress_snps int,
                         track_len int,
                         cum_cov int,
                         total_tracks int,
                         index(start, end),
                         index(ind))""" % (tracks_tmp_filtered_table))             

  tmp_cur.execute("""insert into %s
                     select * from %s
                     where ind = '%s'
                     and arm = '%s'
                     order by ict_id""" % (tracks_tmp_filtered_table, tracks_table, ind, arm))

  tmp_cur.execute("""delete from %s
                     where ind = '%s'
                     and arm = '%s'""" % (tracks_table, ind, arm))
  cnx.commit() 
   
  merge_err_cursor = cnx.cursor()
  merge_err_cursor.execute("""select straight_join t.ict_id + 1, t3.ict_id - 1, 'err'
                              from %s t, %s t3, %s t2
                              where t3.ict_id between t.ict_id + 2 and t.ict_id + 102
                              and t2.ict_id between t.ict_id + 1 and t3.ict_id - 1
                              and t3.genotype not like 'err%%'
                              and t.genotype not like 'err%%'
                              group by t.ict_id, t3.ict_id, t3.genotype
                              having sum(if(t2.genotype not like 'err%%', 1, 0)) = 0
                              """ % (tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))


  filtered_ict_ids = {}
  last_ict_id_start = -99
  last_ict_id_end = -99   
  for (ict_id_start, ict_id_end, genotype) in merge_err_cursor:
    if ict_id_start == last_ict_id_end:
      filtered_ict_ids[last_ict_id_start]['end'] = ict_id_end
    else:
      filtered_ict_ids[ict_id_start] = {'end': -1, 'genotype': 'zzz'}
      filtered_ict_ids[ict_id_start]['end'] = ict_id_end
      filtered_ict_ids[ict_id_start]['genotype'] = genotype
      last_ict_id_start = ict_id_start
      
    last_ict_id_end = ict_id_end
  
  err_filter_tracks('merge_err', tracks_tmp_filtered_table, tracks_table, arm, ind, filtered_ict_ids, hmm_table)
  tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
  cnx.commit()  


   


#begin
(hmm_table, repeats_table) = sys.argv[1:]
tracks_table = hmm_table + '_tracks'
tracks_unfiltered_table = hmm_table + '_unftrk'
filter_log_table = hmm_table + '_fillog'
tmp_repeat_table = hmm_table + '_tmprep'

print(tmp_repeat_table, tracks_table, repeats_table)

create_cur = cnx.cursor()

create_cur.execute("""drop table if exists %s""" % (tracks_unfiltered_table))
create_cur.execute("""create table %s
                      (ict_id int auto_increment primary key,
                       ind varchar(30),
                       genotype varchar(20),
                       arm varchar(5),
                       start int,
                       end int,
                       total_snps int,
                       introgress_snps int,
                       track_len int,
                       cum_cov int,
                       index(start, end),
                       index(ind))""" % (tracks_unfiltered_table))             

create_cur.execute("""drop table if exists %s""" % (tracks_table))
create_cur.execute("""create table %s
                      (ict_id int auto_increment primary key,
                       ind varchar(30),
                       genotype varchar(20),
                       arm varchar(5),
                       start int,
                       end int,
                       total_snps int,
                       introgress_snps int,
                       track_len int,
                       cum_cov int,
                       per_repeat decimal(5,2),
                       index(start, end),
                       index(ind))""" % (tracks_table))             

create_cur.execute("""drop table if exists %s""" % (filter_log_table))
create_cur.execute("""create table %s
                      (fl_id int auto_increment primary key,
                       ind varchar(30),
                       filter_name varchar(30),
                       tracks int,
                       genotype varchar(20),
                       index(filter_name))""" % (filter_log_table))             

create_cur.close()


arm_cursor = cnx.cursor()
arm_cursor.execute("""select distinct arm
                      from %s
                      order by arm""" % (hmm_table))

for (arm,) in arm_cursor:
  print "\t".join([datetime.now().strftime("%y-%m-%d %H:%M:%S"), arm])
  ind_cursor = cnx.cursor()
  ind_cursor.execute("""select distinct ind
                        from %s
                        where arm = '%s'""" % (hmm_table, arm))
  for (ind,) in ind_cursor:
    snp_cursor = cnx.cursor()
    snp_cursor.execute("""select arm, pos, 
                          genotype, pos, a_cov + b_cov
                          from %s
                          where arm = '%s'
                          and ind = '%s'
                          order by pos""" % (hmm_table, arm, ind))
    
    start = 0
    last_pos = 0    
    last_genotype = '' 
    total_snps = 0
    introgress_snps = 0    
    cum_cov = 0             
    for (arm, pos, genotype, cov) in snp_cursor:
      if genotype != last_genotype and last_genotype != '':
        insert_cursor = cnx.cursor()
        insert_cursor.execute("""insert into %s
                                 values
                                 (null, '%s', '%s', '%s', %s, %s, %s, %s, %s, %s)""" % (tracks_unfiltered_table, ind, last_genotype, arm, min(start, last_pos), max(start, last_pos), total_snps, introgress_snps, last_pos - start + 1, cum_cov))
        cnx.commit()                               
        insert_cursor.close()
        start = pos
        total_snps = 0
        introgress_snps = 0
        cum_cov = 0
      
      if start == 0:
        start = pos
        
      last_genotype = genotype
      last_pos = pos
      total_snps += 1
      cum_cov += cov
      if genotype != 'homo_1' and genotype != 'err_homo_1':
        introgress_snps += 1
    
    snp_cursor.close()
    insert_cursor = cnx.cursor()
    insert_cursor.execute("""insert into %s
                             values
                             (null, '%s', '%s', '%s', %s, %s, %s, %s, %s, %s)""" % (tracks_unfiltered_table, ind, last_genotype, arm, min(start, last_pos), max(start, last_pos), total_snps, introgress_snps, last_pos - start + 1, cum_cov))

    insert_cursor.execute("""insert into %s
                             select null, ind, genotype, arm, start, end,
                                    total_snps, introgress_snps, track_len, cum_cov, 0
                             from %s
                             where ind = '%s'
                             and arm = '%s'
                             order by ict_id""" % (tracks_table, tracks_unfiltered_table, ind, arm))
    cnx.commit()                         
    insert_cursor.close()    
  
    tracks_tmp_filtered_table = tracks_unfiltered_table[:-6] + 'tmperr'
    
    tmp_cur = cnx.cursor()
    tmp_cur.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
    tmp_cur.execute("""create table %s
                       (ict_id int primary key auto_increment,
                        ind varchar(30),
                        genotype varchar(20),
                        arm varchar(5),
                        start int,
                        end int,
                        total_snps int,
                        introgress_snps int,
                        track_len int,
                        cum_cov int,
                        total_tracks int,
                        index(start, end),
                        index(ind))""" % (tracks_tmp_filtered_table))             
    
    tmp_cur.execute("""insert into %s
                       select t.*, 1 
                       from %s t
                       where ind = '%s'
                       and arm = '%s'
                       order by ict_id""" % (tracks_tmp_filtered_table, tracks_unfiltered_table, ind, arm))
    
    cnx.commit() 

    
#filter tracks    
    merge_err_filter(tracks_tmp_filtered_table, ind, arm, hmm_table)
    merge_het_homo_2_tmp_filter(tracks_tmp_filtered_table, ind, arm, hmm_table)
    
    err_cursor = cnx.cursor()
    tracks_err_convert_table = tracks_unfiltered_table[:-6] + 'to_err'
    err_cursor.execute("""drop table if exists %s""" % (tracks_err_convert_table))
    err_cursor.execute("""create table %s
                          (te_id int primary key auto_increment,
                           start int,
                           end int,
                           index(start))""" % (tracks_err_convert_table))             
    
    
    err_cursor.execute("""insert into %s
                          select straight_join null, t2.start, t2.end
                          from %s t2, %s t, %s t3
                          where t.ict_id = t2.ict_id - 1
                          and t3.ict_id = t2.ict_id + 1
                          and t2.genotype in ('homo_2', 'het')
                          and t2.total_snps < 15
                          and t2.introgress_snps / (if(t.genotype = 'err', t.total_snps, 0) + if(t3.genotype = 'err', t3.total_snps, 0)) < 3
                          order by t2.start""" % (tracks_err_convert_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table, tracks_tmp_filtered_table))
    
    err_cursor.execute("""update %s t, %s e
                          set t.genotype = concat('err_', t.genotype)
                          where t.start = e.start
                          and arm = '%s'""" % (tracks_table, tracks_err_convert_table, arm))
    
    
    
    err_cursor.execute("""drop table if exists %s""" % (tracks_tmp_filtered_table))
    err_cursor.execute("""drop table if exists %s""" % (tracks_err_convert_table))

    singleton_filter(tracks_table, ind, arm, hmm_table)
    merge_filter(tracks_table, ind, arm, hmm_table)
    
    remove_err_filter(tracks_table, ind, arm, hmm_table, 100)
    #check to see if any error states remain
    #if so, boost the err_runs_filter search size and rerun all filters
    #the larger search size slows down the query by an oder of magnitude, so I only do this as a last resort
    err_check_cursor = cnx.cursor()
    err_check_cursor.execute("""select count(*)
                                from %s
                                where genotype like 'err%%'""" % (tracks_table))
                                
    for (num_err_tracks,) in arm_cursor:
      if num_err_tracks > 0:
        remove_err_filter(tracks_table, ind, arm, hmm_table, 200)

    
    err_check_cursor.execute("""select count(*)
                                from %s
                                where genotype like 'err%%'""" % (tracks_table))
                                
    for (num_err_tracks,) in arm_cursor:
      if num_err_tracks > 0:
        remove_err_filter(tracks_table, ind, arm, hmm_table, 400)
    
    
    merge_small_homo_1_filter(tracks_table, ind, arm, hmm_table)
    merge_het_homo_2_filter(tracks_table, ind, arm, hmm_table)




print(tmp_repeat_table, tracks_table, repeats_table)
#update per_repeat
repeat_cursor = cnx.cursor()
repeat_cursor.execute("""create table %s
                         as select t.ict_id, sum(rep_len) / (end - start + 1) * 100 per_repeat
                         from (select straight_join ict_id, r.end - t.start + 1 rep_len
                               from %s t, %s r
                               where t.arm = r.arm
                               and r.start < t.start
                               and r.end between t.start and t.end
                               union all
                               select straight_join ict_id, t.end - r.start + 1
                               from %s t, %s r
                               where t.arm = r.arm
                               and r.start between t.start and t.end
                               and r.end > t.end
                               union all
                               select straight_join ict_id, r.end - r.start + 1
                               from %s t, %s r
                               where t.arm = r.arm
                               and r.start between t.start and t.end
                               and r.end between t.start and t.end
                               union all
                               select straight_join ict_id, t.end - t.start + 1
                               from %s r, %s t
                               where t.arm = r.arm
                               and t.start between r.start and r.end
                               and t.end between r.start and r.end
                               and t.start > r.start
                               and t.end < r.end) x, %s t
                         where t.ict_id = x.ict_id
                         group by t.ict_id, start, end""" % (tmp_repeat_table, tracks_table, repeats_table, tracks_table, repeats_table, tracks_table, repeats_table, repeats_table, tracks_table, tracks_table))

repeat_cursor.execute("""alter table %s add unique index(ict_id)""" % (tmp_repeat_table))


repeat_cursor.execute("""update %s t, %s x
                         set t.per_repeat = x.per_repeat
                         where t.ict_id = x.ict_id""" % (tracks_table, tmp_repeat_table))

repeat_cursor.execute("""drop table %s""" % (tmp_repeat_table))

repeat_cursor.close()



#check for track filtering problems
problem_cursor = cnx.cursor()
problem_cursor.execute("""select t.ind
                          from (select ind, sum(total_snps) snps
                                from %s 
                                group by ind) u,
                               (select ind, sum(total_snps) snps
                                from %s 
                                group by ind) t
                          where t.ind = u.ind
                          and t.snps <> u.snps
                          order by t.ind""" % (tracks_unfiltered_table, tracks_table))


print("")
no_problems = 'Y'
for (ind) in problem_cursor:
  if no_problems == 'Y':
    print('Problems encountered applying track filtering rules:')
  
  print(ind)
  no_problems = 'N'

if no_problems == 'Y':
  print('No track filtering problems encountered')

print("\n")

