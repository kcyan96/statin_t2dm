library(data.table)
library(readxl)
library(tableone)
library(MatchIt)
library(survival)
library(pbapply)

M <- function(ymd) fifelse(substr(ymd,5,6)=="00", as.integer(NA), as.integer(substr(ymd,1,4))*12 + as.integer(substr(ymd,5,6)))
until <- function(panel, cond) panel[, shift(cummax(eval(parse(text=cond))), fill=0), id][, V1==0]
fmt <- function(x,n) format(round(x,n), nsmall=n)
merge_ <- function(L, by=NULL, all=F, ...) Reduce(function(x,y) merge(x, y, by, all=all, sort=F, ...), L)
as.data.table.glm <- function(m, robust=F, cluster=NULL, vcov.=NULL) {
  if(!is.null(vcov.)) m.vcov <- vcov.
  else if(!is.null(cluster)) m.vcov <- multiwayvcov::cluster.vcov(m, m$data[[cluster]])
  else if(robust) m.vcov <- sandwich::vcovHC(m, type="HC1")
  else m.vcov <- vcov(m)
  t <- lmtest::coeftest(m, vcov.=m.vcov)
  data.table(fmt(cbind(Estimate=exp(t[,"Estimate"]), exp(lmtest::coefci(m, vcov.=m.vcov)), p=t[,"Pr(>|z|)"]), 3), keep.rownames=T)
}

panel_data <- readRDS("../cohort_0516/6.panel_data.RDS")
covar_def <- setDT(read_excel("codes.xlsx", sheet="covariates"))[!is.na(Name)]
matched <- list()
for(MONTH in format(seq(as.Date("2005-01-01"), as.Date("2016-12-31"), by="month"), "%Y%m")) {
  ct <- panel_data[month==MONTH][(M(month)-M(regdate))>=12][dead==0 & xfrd==0][DM==1][age>=25 & age<=84][dx.chd==0 & dx.mi==0 & dx.stroke==0 & dx.hf==0][dx.myopathy==0 & dx.liverdz==0][dx.rhd==0 & dx.cancer==0 & dx.schizo==0 & rx.antipsychotic==0][!grepl("1", substr(rx.statin.map, m_erx-24, m_erx-1))][rx.lipid_other==0 & rx.aspirin==0 & rx.digoxin==0][, complete := as.integer(!is.na(smoking) & smoking!="smoker (unknown quantity)" & !is.na(bmi) & !is.na(sbp) & !is.na(sbp_sd) & !is.na(chol) & !is.na(hdl) & !is.na(ldl) & !is.na(tg) & !is.na(hba1c) & !is.na(egfr))][complete==1]  
  ct[, tc_hdl := chol/hdl][, nhdl := chol-hdl][, qrisk := QRISK3(age, sex, ethnicity, smoking, rep.int(0,nrow(ct)), rep.int(1,nrow(ct)), dx.fh_cvd, dx.ckd_345, dx.af, dx.htn & (rx.ras | rx.bb | rx.ccb | rx.diuretic | rx.htn_other), dx.migraine, dx.ra, dx.sle, dx.smi, rx.antipsychotic2, rx.steroid, dx.erect | rx.erect, chol/hdl, sbp, sbp_sd, bmi, townsend)][, qrisk_cat := cut(qrisk, c(0,10,20,30,100), right=F, include.lowest=T)][, ldl_cat := cut(ldl, c(0,1.8,2.6), right=F, include.lowest=F)][, nhdl_cat := cut(nhdl, c(0,2.6,3.4), right=F, include.lowest=F)][, risk_cat := factor(paste0("QRISK ", qrisk_cat))]
  ct[, statin := as.integer(grepl("1", substr(rx.statin.map, m_erx, m_erx)))]
  set.seed(475); M.glm <- lapply(1:nlevels(ct$risk_cat), function(n) if(ct[risk_cat==levels(risk_cat)[n], uniqueN(statin)]==2) matchit(as.formula(paste0("statin ~ ", paste0(covar_def[Baseline==1, Name], collapse=" + "))), data=ct[risk_cat==levels(risk_cat)[n]], ratio=4, caliper=0.2) else NULL)
  matched[[MONTH]] <- lapply(1:nlevels(ct$risk_cat), function(n) if(!is.null(M.glm[[n]])) ct[risk_cat==levels(risk_cat)[n]][M.glm[[n]]$weights > 0] else NULL)
  for(i in 1:nlevels(ct$risk_cat)) {ct <- matched[[MONTH]][[i]]}; rm(i)
  names(matched[[MONTH]]) <- levels(ct$risk_cat); rm(ct, M.glm)
}; rm(MONTH)
cohort.matched <- setnames(rbindlist(lapply(1:uniqueN(unlist(lapply(matched, names))), function(q) rbindlist(lapply(matched, function(md) md[[sort(unique(unlist(lapply(matched, names))))[q]]]))), idcol="strata"), "month", "index.month")[, id := paste0(strata,"_",index.month,"_",uid)]
tb1 <- lapply(1:max(cohort.matched$strata), function(s) if(cohort.matched[strata==s, .N]>0) table1(cohort.matched[strata==s]) else NULL)

ipw_data <- setorder(merge(cohort.matched[, .(id, uid, index.month)], panel_data, by="uid")[month >= index.month], id, month)[, smoking := fifelse(smoking=="smoker (unknown quantity)", as.integer(NA), as.integer(smoking))][, smoking := nafill(smoking, type="locf"), id][, smoking := factor(smoking, 1:5, labels=c("non-smoker","ex-smoker","light smoker","moderate smoker","heavy smoker"))][, qrisk := QRISK3(age, sex, ethnicity, smoking, rep.int(0,nrow(ipw_data)), rep.int(1,nrow(ipw_data)), dx.fh_cvd, dx.ckd_345, dx.af, dx.htn & (rx.ras | rx.bb | rx.ccb | rx.diuretic | rx.htn_other), dx.migraine, dx.ra, dx.sle, dx.smi, rx.antipsychotic2, rx.steroid, dx.erect | rx.erect, chol/hdl, sbp, sbp_sd, bmi, townsend, override=T)][, statin_trmt := as.integer(!is.na(rx.statin.map) & (substr(rx.statin.map, m_erx, m_erx)=="1" | (function(r) !is.na(r) & r<=12)(m_erx - stringi::stri_locate_last_fixed(substr(rx.statin.map, 1, m_erx-1),"1")[,1] + stringi::stri_locate_first_fixed(substr(rx.statin.map, m_erx+1, nchar(rx.statin.map)),"1")[,1] -1)))][, oc.cvd := as.integer(oc.mi | oc.stroke | oc.hf)][, oc.lfabn := cummax(lab.lfabn), id]
OUTCOMES <- c("dead","oc.cvd","oc.mi","oc.hf","oc.stroke","oc.stroke_isch","oc.stroke_haem","oc.myopathy","oc.myopathy_myalgia","oc.liverdz","oc.lfabn")
ipw_data <- setorder(merge(setnames(cohort.matched[, c("strata", "risk_cat", "id", "index.month", "statin", covar_def[Baseline==1, Name]), with=F], covar_def[Baseline & Varying, Name], paste0(covar_def[Baseline & Varying, Name], ".bl")), ipw_data[, c("id", "month", "xfrd", grep("statin_trmt", names(ipw_data), value=T), unique(c(covar_def[Varying==1, Name], "chol", "lab.lfabn", OUTCOMES))), with=F], by="id"), strata, id, month)[, deviate_any := cummax(as.integer(statin_trmt != statin)), id][until(ipw_data, "deviate_any"), deviate := fifelse(statin==1, as.integer(deviate_any & !oc.myopathy_myalgia & !oc.liverdz & !oc.lfabn), as.integer(deviate_any & qrisk<20))][, deviate := nafill(deviate, type="locf"), id][, t := M(month)-M(index.month)][, t2 := t^2]
gen_ipw <- function(pd, oc, var) {
  ipw_spec <- list(n1 = list("statin==1", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)]),
                   d1 = list("statin==1", c(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], covar_def[Varying==1, Name])),
                   n0 = list("statin==0", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)]),
                   d0 = list("statin==0", c(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], covar_def[Varying==1, Name])))
  ipw_end <- list(deviate = paste0(oc, "|dead|xfrd|deviate_any"), xfrd = paste0(oc, "|dead|xfrd"))[[var]]
  pd[, in_ipw_model := as.integer(until(pd, ipw_end) & t!=0)]  
  ipw_model <- pblapply(ipw_spec, function(s) glm(as.formula(paste0(var, " ~ t + t2 + ", paste0(s[[2]], collapse=" + "))), data = pd[in_ipw_model==1][eval(parse(text=s[[1]]))], family = binomial(link="logit"), model=F, y=F))
  pd[in_ipw_model==1 & statin==1, `:=`(n = ipw_model$n1$fitted.values, d = ipw_model$d1$fitted.values)][in_ipw_model==1 & statin==0, `:=`(n = ipw_model$n0$fitted.values, d = ipw_model$d0$fitted.values)]
  pd[, w := fifelse(in_ipw_model==1, (1-n)/(1-d), 1)]                                             # Fitted probs of censored, but need probs of not being censored, so 1-prob for everyone
  pd[, ws := pmin(pmax(w, quantile(w, 0.01)), quantile(w, 0.99))][, wt := cumprod(ws), id][, wts := pmin(pmax(wt, quantile(wt, 0.01)), quantile(wt, 0.99))]
  return(pd$wts)
}
for(s in sort(unique(ipw_data$strata))) {for(oc in OUTCOMES) {for(typ in c("deviate","xfrd")) {ipw_data[strata==s, c(paste0("wts_", typ, ".", oc)) := gen_ipw(ipw_data[strata==s], oc, typ)]; gc()}; rm(typ)}; rm(oc)}; rm(s)
for(oc in OUTCOMES) ipw_data[, c(paste0("wts.", oc)) := get(paste0("wts_deviate.", oc)) * get(paste0("wts_xfrd.", oc))]; rm(oc)

est_fmla <- paste0(" ~ statin + t + t2 + ", paste0(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], collapse=" + "))
popu_data <- ipw_data[t==0, .SD, .SDcols = c("id", "strata", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)])]
popu_data <- setorder(rbindlist(lapply(0:119, function(n) copy(popu_data)[, t := n]))[, t2 := t^2], id, t)
res_itt <- pblapply(sort(unique(ipw_data$strata)), function(s) pbsapply(OUTCOMES, function(oc) {  
  est_data <- ipw_data[strata==s & until(ipw_data, paste0(oc,"|dead|xfrd")), .SD, .SDcols = c("id", "statin", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts_xfrd.", oc))]
  evt <- est_data[, .(N=uniqueN(id), Events=sum(get(oc)), FU_median=as.numeric(median(.SD[,.N,id]$N)), FU_total=.N), keyby=-statin][, lapply(.SD, function(c) paste0(c, collapse="/"))][, statin := NULL]
  m <- glm(as.formula(paste0(oc, est_fmla)), data = est_data, weights = est_data[, get(paste0("wts_xfrd.", oc))], family = quasibinomial(link="logit"))
  res <- setDT(cbind(evt, as.data.table(m, robust=T)[rn=="statin"]))
print(res); rm(m, est_data); gc(); return(res)}, simplify=F))
res_itt <- rbindlist(pblapply(res_itt, function(x) rbindlist(x, id="Outcome")[, rn := NULL]), id="strata")[, strata := levels(ipw_data$risk_cat)[strata]]
abs_risk_itt <- pblapply(sort(unique(ipw_data$strata)), function(s) pbsapply(OUTCOMES, function(oc) {cat("\n", s, oc, "...\n")
  est_data <- ipw_data[strata==s & until(ipw_data, paste0(oc,"|dead|xfrd")), .SD, .SDcols = c("id", "statin", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts_xfrd.", oc))]
  m <- glm(as.formula(paste0(oc, paste0(est_fmla, " + I(statin * t) + I(statin * t2)"))), data = est_data, weights = est_data[, get(paste0("wts_xfrd.", oc))], family = quasibinomial(link="logit"))
  arc <- merge_(lapply(c(1,0), function(treat) setnames(cbind(popu_data[strata==s, .(id, t)], prob = predict(m, copy(popu_data[strata==s])[, statin := treat], type = "response"))[, surv := cumprod(1-prob), id][, mean(1-surv), t], "V1", if(treat==1) "treated" else "control")), by="t")
rm(m, est_data); gc(); return(arc)}, simplify=F))
abs_risk_itt <- melt(rbindlist(pblapply(abs_risk_itt, function(x) rbindlist(x, id="Outcome")), id="strata")[, `:=`(diff = treated-control, ratio = treated/control)], id.vars=c("strata","Outcome","t"), variable.name="Type", value.name="Estimate")[, strata := levels(ipw_data$risk_cat)[strata]]
res_pp <- pblapply(sort(unique(ipw_data$strata)), function(s) pbsapply(OUTCOMES, function(oc) {cat("\n", s, oc, "...\n")
  est_data <- ipw_data[strata==s & until(ipw_data, paste0(oc,"|dead|xfrd|deviate")), .SD, .SDcols = c("id", "statin", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts.",oc))]
  evt <- est_data[, .(N=uniqueN(id), Events=sum(get(oc)), FU_median=as.numeric(median(.SD[,.N,id]$N)), FU_total=.N), keyby=-statin][, lapply(.SD, function(c) paste0(c, collapse="/"))][, statin := NULL]
  m <- glm(as.formula(paste0(oc, est_fmla)), data = est_data, weights = est_data[, get(paste0("wts.", oc))], family = quasibinomial(link="logit"))
  res <- setDT(cbind(evt, as.data.table(m, robust=T)[rn=="statin"]))
print(res); rm(m, est_data); gc(); return(res)}, simplify=F))
res_pp <- rbindlist(pblapply(res_pp, function(x) rbindlist(x, id="Outcome")[, rn := NULL]), id="strata")[, strata := levels(ipw_data$risk_cat)[strata]]
abs_risk_pp <- pblapply(sort(unique(ipw_data$strata)), function(s) pbsapply(OUTCOMES, function(oc) {cat("\n", s, oc, "...\n")
  est_data <- ipw_data[strata==s & until(ipw_data, paste0(oc,"|dead|xfrd|deviate")), .SD, .SDcols = c("id", "statin", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts.",oc))]
  m <- glm(as.formula(paste0(oc, paste0(est_fmla, " + I(statin * t) + I(statin * t2)"))), data = est_data, weights = est_data[, get(paste0("wts.", oc))], family = quasibinomial(link="logit"))
  arc <- merge_(lapply(c(1,0), function(treat) setnames(cbind(popu_data[strata==s, .(id, t)], prob = predict(m, copy(popu_data[strata==s])[, statin := treat], type = "response"))[, surv := cumprod(1-prob), id][, mean(1-surv), t], "V1", if(treat==1) "treated" else "control")), by="t")
rm(m, est_data); gc(); return(arc)}, simplify=F))
abs_risk_pp <- melt(rbindlist(pblapply(abs_risk_pp, function(x) rbindlist(x, id="Outcome")), id="strata")[, `:=`(diff = treated-control, ratio = treated/control)], id.vars=c("strata","Outcome","t"), variable.name="Type", value.name="Estimate")[, strata := levels(ipw_data$risk_cat)[strata]]
