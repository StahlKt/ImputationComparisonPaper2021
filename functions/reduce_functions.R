#reduce_functions

wrap_reduceResults<-function(id.list,reg.result){
  if(length(id.list)==1){
    build_reduced_tables(reduceResults(fun=rbind,
                                       init=id.list[1],
                                       reg=reg.result)) 
  }
  else{
   build_reduced_tables(reduceResults(fun=rbind,
                                     ids=id.list,
                                     reg=reg.result))
  }

}

wrap_reduceResults_pbwt<-function(id){
  build_reduced_tables(reduceResults(fun = rbind,
                                     init = id,
                                     reg = tmp_results))
}

get_job_id_only<-function(in.table){
  in.table$job.id
}



merge_reduced_tables<-function(list1,list2){
  
    list_tables<-list(rbind(list1[[1]],list2[[1]]),
                      rbind(list1[[2]],list2[[2]]),
                      rbind(list1[[3]],list2[[3]]))
          
}
