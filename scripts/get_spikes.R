# faire sur le flotteurs de l'Atlantique Nord avec du zoo
# 6901180
# 6901181


get_spikes <- function(x, p, xerr=NA, pr=NA, max_iter=NA){
  require(pracma)
  require(arules)
# %SPIKE_ANALYSIS a point is considered as a spike if the difference between x and a 15 points Hampel
# %       filter of x is greater than 3 times the scaled Median Absolute Deviation
# %       also return negative spikes separately
# %
# % INPUTS (units do not matter):
#   %     x <Nx1 double> profile of bio-optical data, for example:
#   %         - particulate backscattering (beta, bbp)
# %         - fluorescent dissolved organic matter (fdom)
# %         - chlorophyll a fluorescence  (fchl)
# % OPTIONAL INPUT:
#   %     p <Nx1 double> pressure of profile (add sampling frequency check)
# %     xerr <1x1 double> absolute threshold for spike detection
# %         default: 3x the scaled Median Absolute Deviation (MAD) of entire profile (non-parametric)
# %         recommended: for fdom: ???; for bbp: ???; for fchl: ???;                 (parametric)
# %     pr <1x1 double> minimum sampling resolution threshold
# %         default: 3x the median sampling resolution of entire profile             (non-parametric)
# %         recommended: 5 dBar (typically no spike layers are found below 2.25 dBar)(parametric)
# %     max_iter <1x1 integer> maximum number of intereation to run Hampel Filter
# %         set to 0 to run the Hampel filter only once
# %         default: 100
# %     
# %
# % OUTPUTS:
#   %     spikes <NX1 boolean> return true for every indices containing a positive spike
# %     neg_spikes <NX1 boolean> return true for every indices containing a negative spike
# %     spikes_p <NX1 boolean> pressure of each positive spikes
# %     neg_spikes_p <NX1 boolean> pressure of each negative spikes
# %
# % References:
#   %   Rousseeuw, P. J., and C. Croux (1993), Alternatives to the Median Absolute Deviation,
# %     Journal of the American Statistical Association, 88(424), 1273?1283.
# %
# % author: Nils Haentjens
# % created: May 1, 2018
  x <- tibble(prof = p, x = x) %>% arrange(prof) %>% pull(x)
  p <- tibble(prof = p, x = x) %>% arrange(prof) %>% pull(prof)

    DEBUG = F

# Check input
    #nargin <- length(as.list(match.call())) -1
    nargin = 2
    
    if (nargin < 2 | length(p) == 0){p_flag = F
    }else {p_flag = T}
    
    if (nargin < 3 | length(xerr) == 0){ xerr_flag = F
                                         xerr <- NA
    }else {xerr_flag = T}
    
    if (nargin < 4 | length(pr) == 0){ pr_flag = F
    }else {p_flag = T} 
      #ignored if ~p_flag
    
    if (nargin < 5 | length(max_iter) == 0){max_iter = 100}
    
    
    if (!is.vector(x)) {warning('x is not a row vector (Nx1).')}
    if (p_flag & any(length(p) != length(x))) {warning('Input vectors are different size.')}

# # % Check output
#     if nargout > 2 || (nargout == 2 && ~p_flag); neg_spike_flag = true; else; neg_spike_flag = false; end
# 
# # % Live debug
#     if (DEBUG)
#     fig(1);
#     subplot(1,3,2); hold('on');
#     plot(x,p,'.-'); 
#     end

# % Prepare input
    spikes <- logical(length = length(x))  #attention il faut que toutes les valeurs soient = à FALSE

    if (p_flag==T){
    # ignore nan of x and p for processing
          sel = which(!is.na(x) & !is.na(p) & !is.infinite(x) & !is.infinite(p))
    
      
          if (length(sel) < 3) {
            warning('GET_SPIKES:EMPTY_INPUT', 'get_spikes: Not enough valid values.')
          
            if (nargout == 2 & !p_flag) {varargout = spikes}
            else {varargout <- list(vector(), spikes, vector())}
            return(varagout)
          }
          
          p1 = p[sel]
          # ignore low sampling resolution of profile
          delta <- matrix(nrow=length(p1), ncol=1)
          delta[1,1] <-  abs(diff(p1[1:2]))
          delta[2:(length(p1)-1),1] <- abs((p1[3:length(p1)] - p1[1:(length(p1)-2)]) / 2)
          delta[length(p1),1] <- abs(diff(p1[(length(p1)-1):length(p1)]))  # / dBar 
          if (!pr_flag){ 
              pr = 3 * median(delta)
              subsel = delta < pr
              sel = sel[subsel]
          }
          if (DEBUG){
          #subplot(1,3,1)
              split.screen(c(1,3)) 
              plot(delta, p)
              set(gca, 'ydir', 'reverse') #renverser l'axe Y
          }
          # % keep only relevant sampling
          x = x[sel]
          if(length(p)!=length(sel)){
            p = p[sel]
          }else{
            p=p1
          }
      }else{
        # % only x
          sel = which(!is.na(x) & !is.infinite(x))
          x = x[sel]
        }
  
    if(length(x)>=(35*2)){
        # % Apply Hampel filter (same as median filter except onlx change values of spikes)
        xd <- hampel(x, 35) # > 5 is good
        # % Recursive Hampel filter
        xdo <- x
        i <- 0
        while (any(xd$y != xdo) & i < max_iter){
            xdo <- xd$y
            xd <- hampel(xd$y, 35)
            i <- i + 1
        }
        if (i == max_iter){warning('GET_SPIKES:MAX_ITER', 'get_spikes: Maximum iteration (%d) reached.', max_iter)}
        
        if (!xerr_flag){
            # % Compute Scaled Median Absolute Deviation (MAD)
            smad = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(x-median(x)))
            # Noisy signal (no significant spike)
            if (smad == 0){
            # % format output & return
              if (nargout == 2 & !p_flag){varargout = {spikes}}
              else {varargout = list(vector(), spikes, vector()) } 
            return(varargout)
            }
    # % Set spike detection threshold to 3 scaled MAD
        xerr = 3 * smad
        }
    
    # % Get spikes
        x_delta <-  x - xd$y
        spikes[sel] = x_delta > xerr
    
    # % Live debug
    #     if (DEBUG){
    #     subplot(1,3,2)
    #     plot(xd, p, '.-')
    #     plot(x(spikes(sel)), p(spikes(sel)), 'o')
    #     set(gca, 'ydir', 'reverse')
    #     y_lim = ylim()
    #     subplot(1,3,3); hold('on')
    #     plot(x_delta, p, '.-')
    #     plot(xerr .* ones(2,1), ylim(), 'k', 'LineWidth', 2)
    #     ylim(y_lim)
    #     set(gca, 'ydir', 'reverse')
    #     xlim([0 0.01])
    # # %   set(gca, 'xscale', 'log');
    #     drawnow()
    #     waitforbuttonpress()
    # # %   pause(1);
    #     }
    # 
    # # % Get negative spikes
    #     if (neg_spike_flag){
    #     neg_spikes = false(size(spikes))
    #     neg_spikes(sel) = -x_delta > xerr
    #     }
    # 
    # % Set output
        if (p_flag){
          spikes <- which(spikes[sel]==T)
          
        }
        
        if(length(spikes!=0)){
        x_spikes <- x[spikes] # profil déspiké
        p_spikes <- p[spikes] # profondeurs ou il y a les spikes
        }
        return(list(p_spikes=p_spikes, x_spikes=x_spikes ))
    
    }else{
      return(list(p_spikes="sample is too small", x_spikes="sample is too small"))
    }
}


file <- read_csv("data/qcCSV/6901180.csv")
profil <- file[which(file$cycleNumber==23),]

x <- profil$BBP700
p <- profil$PRES
p_spikes <- get_spikes (x,p)

# x1 <- profil$BBP700_ORIG
# x <- x1

plot(x,-p)
plot(x2)





