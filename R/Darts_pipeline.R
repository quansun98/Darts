#' The Darts Bayesian Hypothesis Testing function without replicate
#' 
#' It calculates posterior probability of differential splicing, taking input files from rMATs output (or similar file formats). It also provides an option to use iterative empirical Bayesian approach to update the original prior. 
#' @param in_fn Input file name. Required. The input file is in rMATs output format, requiring columns ID,I1,S1,I2,S2,inc_len,skp_len. Additional columns are allowed.
#' @param out_fn Output file name. Required. This function will directly write results to the file. 
#' @param iter Use iteration or not. Logical, TRUE or FALSE. Default TRUE if no informative prior provided. If TRUE, it will start from uninformative prior (if no informative prior provided), and update the prior based on the posterior from the last round, until convergence. Note that if informative prior is provided, it will be set to FALSE with a WARNING message, even if the user specificed TRUE.
#' @param C The cutoff to define significant difference, i.e., delta psi values, where psi is the exon inclusion rate. Default 0.05.
#' @param rescale_meth Rescale method. It will only be considered if informative prior file is provided, taking value 1 or 0. If 1, it will rescale by binary prediction, and by bias otherwise. Default 1.
#' @param rho_fn Informative prior file name. The file format is from Darts_DNN output. If not provided, it will start from uninformative prior. We suggest setting iter=TRUE if no informative prior.
#' @param verbose Default 1.
#' @param epsilon Convergence criteria. It will only on considered if iter=TRUE. Default 0.01.
#' @param random_state random seed. Default 777.
#' @param thread Number of thread to perform parallel computing. Default will use all available threads if not provided.
#'
#' @return A dataframe of posterior probability. The following columns will be added appending to the original input file.
#' \item{rho}{The final prior P1, proportion of differential splicing events, from the last iteration to generate this output.}
#' \item{psi1}{Exon inclusion rate (psi values) for the first group.}
#' \item{psi2}{Exon inclusion rate (psi values) for the second group.}
#' \item{mu.mle}{Maximum likelihood estimator of mean psi value for group 1.}
#' \item{delta.mle}{Maximum likelihood estimator of mean difference in psi values in group 2 compared to group 1.}
#' \item{post_pr}{Posterior probability of this event being differential splicing. We recommend using post_pr > 0.9 to define significant differential splicing, post_pr < 0.2 to define no differential splicing, and in between as unsure.}
#' \item{covg}{Mean coverage for this event.}
#' 
#' @examples 
#' 
#' out = Darts(in_fn = "test_rmats.txt", out_fn = "test_darts.txt", iter=T)
#'
#' @export Darts
#' 

Darts = function(in_fn, out_fn, iter=TRUE, C=0.05, rescale_meth=1, rho_fn=NA, verbose=1, epsilon=0.01, random_state=777, thread=NULL)
{

	if((!is.na(rho_fn)) && iter){
                message("WARNING: informative priors provided, but iteration mode is on, will NOT do iterations")
                iter = FALSE
        }

	if ( is.null(thread ) ) { thread = parallel::detectCores() }
	cat("no. of threads using =", thread, '\n')
	suppressMessages(library('doSNOW'))
	cl = makeSOCKcluster(thread)
	registerDoSNOW(cl)

	#out_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'darts_bht.flat.RData', 'dart_bht.info.RData'))
	#out_table_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'darts_bht.flat.txt', 'dart_bht.info.txt'))
	if(class(in_fn)=="character") {
		data=read.table(in_fn, sep='\t', header=T)
	} else if(class(in_fn)=="data.frame")
	{
		data = in_fn
	}
	if("IJC_SAMPLE_1" %in% colnames(data)) {
		data = match_input_column_names(data)
	}
	res = data
	if(class(rho_fn)=="character")
	{
		rho.df = read.table(rho_fn, header=T)
		if(rescale_meth==1) {
			rho.df$rho = rescale_binary_prediction(rho.df$Y_pred, mu0.new=0.05, mu1.new=0.95, sd0.new=0.1, sd1.new=0.1, boundary=0.05)
		}else
		{
			rho.df$rho = rescale_by_bias(rho.df$Y_pred, rho.df$Y_true)
		}
		idx = match(data$ID, rho.df$ID)
		data$rho = rho.df$rho[idx]
	} else {
		if(! "rho" %in% colnames(data))
			data$rho = rep(0.5, nrow(data))
	}
	#res$mu.mle = NA
	#res$delta.mle = NA
	#res$post_pr = NA
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	delta_quantiles = as.vector(delta_quantiles)
	res_dict = matrix(NA, nrow=nrow(data), ncol=length(delta_quantiles))
	rownames(res_dict) = data$exon_id
	
	pb = utils::txtProgressBar(min = 0, max = nrow(data), initial = 0, style=3)
	progress = function(n) utils::setTxtProgressBar(pb, n)
	opts = list(progress=progress)
	#for(i in 1:nrow(data))

	i_iter = 1
	new_rho = NA
	old_rho = NA

	while(TRUE){

	mp.res = foreach(i=1:nrow(data), .packages='Darts', .options.snow=opts, .combine='rbind') %dopar%
	{
		if(verbose==2 && ! i%%200) cat(format(Sys.time(), "%c"), i,'/',nrow(data), '\n')
		if(verbose==1) utils::setTxtProgressBar(pb,i)
		id=data$ID[i]
		I1=data$I1[i]
		S1=data$S1[i]
		I2=data$I2[i]
		S2=data$S2[i]

		if(!is.na(rho_fn)){
                        this_rho = data$rho[i]
                } else{
                        this_rho = ifelse(is.na(new_rho), 0.5, new_rho)
                }

		inc_len = data$inc_len[i]
		skp_len = data$skp_len[i]
		this_data = list(I1=I1, S1=S1, I2=I2, S2=S2)
		is_healthy = check_data_sanity(this_data)
		if(! is_healthy) {
			this = NULL
		} else {
		set.seed(random_state)
		post = darts_bayesian_gibbs(rho=this_rho, data=this_data, 
			tau1=0.3, tau2=0.03,
			inc_eff_len=inc_len, skp_eff_len=skp_len)
		set.seed(NULL)
		post_cdf = sapply(delta_quantiles, function(x) mean(post[,2]<x))
		res_dict[i,] = post_cdf
		right = which(abs(delta_quantiles-abs(C))<0.0001)
		left = which(abs(delta_quantiles+abs(C))<0.0001)
		
		#res$post_pr[i] = 1-post_cdf[right] + post_cdf[left]
		#res$mu.mle[i] = I1/inc_len / (I1/inc_len + S1/skp_len)
		#res$delta.mle[i] = I2/inc_len / (I2/inc_len + S2/skp_len) - I1/inc_len / (I1/inc_len + S1/skp_len)
		this = data.frame( 
			ID = id, 
			psi1 = round( I1/inc_len / (I1/inc_len + S1/skp_len), 3),
			psi2 = round( I2/inc_len / (I2/inc_len + S2/skp_len), 3),
			delta.mle = round( I2/inc_len / (I2/inc_len + S2/skp_len) - I1/inc_len / (I1/inc_len + S1/skp_len), 4),
			post_pr = round(1-post_cdf[right] + post_cdf[left], 4),
			covg = (I1 + I2 + S1 + S2)/2,
			stringsAsFactors=F
			)
		}
		this
	}

	if(!iter){ 
		message("no iteration performed, finished!")
		break }

	if(is.na(old_rho)){old_rho = 0.5}
	covg_quartile = dplyr::ntile(mp.res$covg, 4) 
    	mp.res_sub = mp.res[which(covg_quartile >= 2),]
    	new_rho = sum(mp.res_sub$post_pr > 0.9) / nrow(mp.res_sub)

        if(new_rho < 0.01){ new_rho = 0.01 
        } else if(new_rho > 0.99) {new_rho = 0.99}

        if(abs(new_rho - old_rho) < epsilon) {
                message(paste("iteration",i_iter,"finished, new prior1 =", new_rho, ", epsilon =", epsilon, ", iteration finished!"))
                break
        } else{
                message(paste("iteration",i_iter,"finished, new prior1 =", new_rho, ", epsilon =", epsilon, ", continue to next iteration..."))
        	i_iter = i_iter + 1
		old_rho = new_rho
        }

        }

	close(pb)
	res = merge(res, mp.res, by='ID')
	#if(!is.null(out_RData_fn)) save(res_dict, file=out_RData_fn)
	write.table(res, file=out_fn, quote=F, row.names=F, sep='\t')
	stopCluster(cl)
	return(res)
}




#' The Darts Bayesian Hypothesis Testing function with replicate
#' 
#' It calculates posterior probability of differential splicing, taking input files from rMATs output (or similar file formats). It also provides an option to use iterative empirical Bayesian approach to update the original prior. 
#' @param in_fn Input file name. Required. The input file is in rMATs output format, requiring columns ID,I1,S1,I2,S2,inc_len,skp_len. Additional columns are allowed. Multiple replicates should be separated by comma within columns I1,S1,I2,S2.
#' @param out_fn Output file name. Required. This function will directly write results to the file. 
#' @param iter Use iteration or not. Logical, TRUE or FALSE. Default TRUE if no informative prior provided. If TRUE, it will start from uninformative prior (if no informative prior provided), and update the prior based on the posterior from the last round, until convergence. Note that if informative prior is provided, it will be set to FALSE with a WARNING message, even if the user specificed TRUE.
#' @param C The cutoff to define significant difference, i.e., delta psi values, where psi is the exon inclusion rate. Default 0.05.
#' @param rescale_meth Rescale method. It will only be considered if informative prior file is provided, taking value 1 or 0. If 1, it will rescale by binary prediction, and by bias otherwise. Default 1.
#' @param rho_fn Informative prior file name. The file format is from Darts_DNN output. If not provided, it will start from uninformative prior. We suggest setting iter=TRUE if no informative prior.
#' @param estim_groupVar_prior If output plot for estimated group variance prior. Default FALSE. If TRUE, if will output a plot with name groupVar.pdf in the working directory.
#' @param is_paired If data are paired between groups. Logical, TRUE or FALSE. Default FALSE.
#' @param pooling If pooling different replicates together as a single input for inference. Logical, TRUE or FALSE. Default FALSE. Suggest setting to FALSE unless necessary.
#' @param verbose Default 1.
#' @param epsilon Convergence criteria. It will only on considered if iter=TRUE. Default 0.01.
#' @param random_state random seed. Default 777.
#' @param thread Number of thread to perform parallel computing. Default will use all available threads if not provided.
#'
#' @return A dataframe of posterior probability. The following columns will be added appending to the original input file.
#' \item{rho}{The final prior P1, proportion of differential splicing events, from the last iteration to generate this output.}
#' \item{psi1}{Exon inclusion rate (psi values) for the first group. Multiple replicates will be separated by comma.}
#' \item{psi2}{Exon inclusion rate (psi values) for the second group. Multiple replicates will be separated by comma.}
#' \item{mu.mle}{Maximum likelihood estimator of mean psi value for group 1.}
#' \item{delta.mle}{Maximum likelihood estimator of mean difference in psi values in group 2 compared to group 1.}
#' \item{post_pr}{Posterior probability of this event being differential splicing. We recommend using post_pr > 0.9 to define significant differential splicing, post_pr < 0.2 to define no differential splicing, and in between as unsure.}
#' \item{covg}{Mean coverage for this event.}
#' 
#' @examples 
#' 
#' out = Darts_replicate(in_fn = "test_rmats_rep.txt", out_fn = "test_darts_rep.txt", iter=T)
#'
#' @export Darts_replicate

Darts_replicate = function(in_fn, out_fn, iter=TRUE, rescale_meth=1, C=0.05, rho_fn=NA, 
	estim_groupVar_prior=FALSE, is_paired=FALSE, pooling=FALSE, epsilon=0.01, verbose=1, random_state=777, thread=NULL)
{

	if((!is.na(rho_fn)) && iter){
                message("WARNING: informative priors provided, but iteration mode is on, will NOT do iterations")
                iter = FALSE
        }

	if ( is.null(thread ) ) { thread = parallel::detectCores() }
	cat("no. of threads using =", thread, '\n')
	suppressMessages(library('doSNOW'))
	cl = makeSOCKcluster(thread)
	registerDoSNOW(cl)

	#out_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'dart_bht.flat.RData', 'darts_bht.info.RData'))
	#out_table_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'dart_bht.flat.txt', 'dart_bht.info.txt'))
	if(class(in_fn)=="character") {
		data=read.table(in_fn, sep='\t', header=T, stringsAsFactors=F)
	} else if(class(in_fn)=="data.frame")
	{
		data = in_fn
	}
	if("IJC_SAMPLE_1" %in% colnames(data)) {
		data = match_input_column_names(data)
	}
	res = data
	if(class(rho_fn)=="character")
	{
		rho.df = read.table(rho_fn, header=T)
		if(rescale_meth==1) {
			rho.df$rho = rescale_binary_prediction(rho.df$Y_pred, mu0.new=0.05, mu1.new=0.95, sd0.new=0.1, sd1.new=0.1, boundary=0.05)
		} else {
			rho.df$rho = rescale_by_bias(rho.df$Y_pred, rho.df$Y_true)
		}
		idx = match(data$ID, rho.df$ID)
		data$rho = rho.df$rho[idx]
	} else {
		if(! "rho" %in% colnames(data))
			data$rho = rep(0.5, nrow(data))
	}
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	delta_quantiles = as.vector(delta_quantiles)
	res_dict = matrix(NA, nrow=nrow(data), ncol=length(delta_quantiles))
	rownames(res_dict) = data[,1]
	
	if(estim_groupVar_prior) {
		prior_fit = estim_group_var.prior_fit(data)
		sgm_prior_fit = prior_fit$fit
		ggplot2::ggsave(file.path('groupVar.pdf'), plot=prior_fit$p, width=7, height=7)
	} else {
			sgm_prior_fit = NA
	}
	
	pb = utils::txtProgressBar(min = 0, max = nrow(data), initial = 0, style=3)
	progress = function(n) utils::setTxtProgressBar(pb, n)
	opts = list(progress=progress)
	#for(i in 1:nrow(data))

	i_iter = 1
	new_rho = NA
	old_rho = NA
	while(TRUE){

	mp.res = foreach(i=1:nrow(data), .packages='Darts', .options.snow=opts, .combine='rbind') %dopar%
	{
		if(verbose==2 && ! i%%200) cat(format(Sys.time(), "%c"), i,'/',nrow(data), '\n')
		if(verbose==1) utils::setTxtProgressBar(pb,i)
		id=data$ID[i]
		I1=data$I1[i]; I1 = as.numeric(strsplit(as.character(I1),',')[[1]])
		S1=data$S1[i]; S1 = as.numeric(strsplit(as.character(S1),',')[[1]])
		I2=data$I2[i]; I2 = as.numeric(strsplit(as.character(I2),',')[[1]])
		S2=data$S2[i]; S2 = as.numeric(strsplit(as.character(S2),',')[[1]])
		
		if(!is.na(rho_fn)){
                        this_rho = data$rho[i]
                } else{
                        this_rho = ifelse(is.na(new_rho), 0.5, new_rho)
                }

		inc_len = data$inc_len[i]
		skp_len = data$skp_len[i]
		
		inc = c(I1, I2)
		skp = c(S1, S2)
		if(is_paired) {
			dat.k=rep(1:length(I1), 2)
		} else {
			dat.k=seq(1,length(inc))
		}
		dat.j = c(rep(1, length(I1)), rep(2, length(I2)))
		this_data = data.frame(k=dat.k, j=dat.j, Inc=inc, Skp=skp)
		set.seed(random_state)
		if(pooling) { 
			this_data = collapse_replicates(this_data)
			is_healthy = check_data_sanity(this_data)
			if(! is_healthy) {
				this = NULL
			} else {
				post = darts_bayesian_gibbs(rho=this_rho, data=this_data, 
					tau1=0.3, tau2=0.03,
					inc_eff_len=inc_len, skp_eff_len=skp_len)
				# res$I1[i] = this_data$I1
				# res$I2[i] = this_data$I2
				# res$S1[i] = this_data$S1
				# res$S2[i] = this_data$S2
			}
		} else {
			is_healthy = check_data_sanity_replicate(this_data)
			if(! is_healthy) {
				this = NULL
			} else {
			post = rdarts_bayesian_gibbs(rho=this_rho, data=this_data, 
				tau1=0.3, tau2=0.03, sgm_prior_fit=sgm_prior_fit,
				inc_eff_len=inc_len, skp_eff_len=skp_len)
			}
		}
		set.seed(NULL)
		if(is_healthy) {
			post_cdf = sapply(delta_quantiles, function(x) mean(post[,2]<x))
			#res_dict[i,] = post_cdf
			right = which(abs(delta_quantiles-abs(C))<0.0001)
			left = which(abs(delta_quantiles+abs(C))<0.0001)
			#res$post_pr[i] = 1-post_cdf[right] + post_cdf[left]
			#res$mu.mle[i] = sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len)
			#res$delta.mle[i] = sum(I2)/inc_len / (sum(I2)/inc_len + sum(S2)/skp_len) - sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len)
			if(pooling) {
				this = data.frame(
					ID = id,
					I1.p = this_data$I1,
					S1.p = this_data$S1,
					I2.p = this_data$I2,
					S2.p = this_data$S2,
					rho = this_rho,
					psi1 = round(sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len), 3),
					psi2 = round(sum(I2)/inc_len / (sum(I2)/inc_len + sum(S2)/skp_len), 3),
					delta.mle = round( sum(I2)/inc_len / (sum(I2)/inc_len + sum(S2)/skp_len) - sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len), 4),
					post_pr = round(1-post_cdf[right] + post_cdf[left], 4),
					covg = (sum(I1) + sum(I2) + sum(S1) + sum(S2))/2,
					stringsAsFactors=F
					)
			} else {
				psi1 = paste( round(I1/inc_len / (I1/inc_len + S1/skp_len), 3), collapse=',')
				psi2 = paste( round(I2/inc_len / (I2/inc_len + S2/skp_len), 3), collapse=',')
				if(is_paired) {
					mu.mle = mean( I1/inc_len / (I1/inc_len + S1/skp_len) )
					delta.mle = mean( I2/inc_len / (I2/inc_len + S2/skp_len) ) - mean( I1/inc_len / (I1/inc_len + S1/skp_len) )
				} else {
					mu.mle = mean( I1/inc_len / (I1/inc_len + S1/skp_len) )
					delta.mle = mean( I2/inc_len / (I2/inc_len + S2/skp_len) - I1/inc_len / (I1/inc_len + S1/skp_len) )
				}
				this = data.frame(
					ID = id,
					rho = this_rho,
					psi1 = psi1,
					psi2 = psi2,
					mu.mle = round(mu.mle, 4),
					delta.mle = round(delta.mle, 4),
					post_pr = round(1-post_cdf[right] + post_cdf[left], 4),
					covg = (mean(I1) + mean(I2) + mean(S1) + mean(S2))/2,
					stringsAsFactors=F
					)
			}
		}
		this
	}

	if(!iter){
                message("no iteration performed, finished!")
                break }

	if(is.na(old_rho)){old_rho = 0.5}

	covg_quartile = dplyr::ntile(mp.res$covg, 4) 
	mp.res_sub = mp.res[which(covg_quartile >= 2),]
    	new_rho = sum(mp.res_sub$post_pr > 0.9) / nrow(mp.res_sub)

        if(new_rho < 0.01){ new_rho = 0.01 
        } else if(new_rho > 0.99) {new_rho = 0.99}

        if(abs(new_rho - old_rho) < epsilon) {
                message(paste("iteration",i_iter,"finished, new prior1 =", new_rho, ", epsilon =", epsilon, ", iteration finished!"))
                break
        } else{
                message(paste("iteration",i_iter,"finished, new prior1 =", new_rho, ", epsilon =", epsilon, ", continue to next iteration..."))
        	i_iter = i_iter + 1
		old_rho = new_rho
        }

        }

	close(pb)
	#if(!is.null(out_RData_fn)) save(res_dict, file=out_RData_fn)
	res = merge(res, mp.res, by='ID')
	stopCluster(cl)
	write.table(res, file=out_fn, quote=F, row.names=F, sep='\t')
	return(res)
}


match_input_column_names = function(df)
{
	rmats_names = c('ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen')
	darts_names = c('ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len')
	name_matching = data.frame(rmats=rmats_names, darts=darts_names)
	idx = match(name_matching$rmats, colnames(df), nomatch=0)
	df = df[,idx]
	idx2 = match(colnames(df), name_matching$rmats)
	colnames(df) = name_matching$darts[idx2]
	return(df)
}


check_data_sanity = function(data)
{
	is_healthy=1
	if(data$I1==0 && data$S1==0) is_healthy=0
	if(data$I2==0 && data$S2==0) is_healthy=0
	#if( min(data$I1, data$S1)<=2 ) is_healthy=0
	#if( min(data$I2, data$S2)<=2 ) is_healthy=0
	#if( sum(data$I1, data$S1)<=20 ) is_healthy=0
	#if( sum(data$I2, data$S2)<=20 ) is_healthy=0
	#if(data$I1==0 && data$I2==0) is_healthy=0
	#if(data$S1==0 && data$S2==0) is_healthy=0
	return(is_healthy)
}

check_data_sanity_replicate = function(data)
{
	is_healthy=1
	n_rep = nrow(data)/2
	if(sum(data$Inc==0 & data$Skp==0)>0) is_healthy=0
	if(all(data$Skp==0)) is_healthy=0
	if(all(data$Inc==0)) is_healthy=0
	#if(all(data$Skp[data$j==1]==0)) is_healthy=0
	#if(all(data$Skp[data$j==2]==0)) is_healthy=0
	#if(all(data$Inc[data$j==1]==0)) is_healthy=0
	#if(all(data$Inc[data$j==2]==0)) is_healthy=0
	return(is_healthy)
}

change_cutoff = function(in_fn, rdata_fn, C)
{
	data=read.table(in_fn, sep='\t', header=T)
	load(rdata_fn)
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	right = which(abs(delta_quantiles-abs(C))<0.0001)
	left = which(abs(delta_quantiles+abs(C))<0.0001)
	post_pr = 1-res_dict[,right] + res_dict[,left]
	data$post_pr = post_pr
	return(data)
}


