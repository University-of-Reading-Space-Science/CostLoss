%a script to produce and validate probabilistic forecasts based on
%different knowledge scenarios. e.g., if we could predict when an ICME
%would arrive, but not its properties; or when CIRs will arrive, etc.
%
%This is a rehash of the python code I wrote a few years back.
%
%(Mathew Owens, 14/11/18)
clear
smjd=date2mjd(1995,1,1);
fmjd=date2mjd(2020,1,1);
percent=0.99; %the percentile to consider

ncat=10; %number of categories for V, B. 4 = quartiles
nbins=10000; %number of bins for CDF [10000]
ncl=100; %number of cost/loss bins
fsize=14; %font size for legends

load([getenv('DBOX'),'matlabdatafiles\omni2_complete_Feb2020.mat']);
load([getenv('DBOX'),'matlabdatafiles\ICMEs_updateSept2019.mat']);

omni=omni_1hour; clear omni_1hour;
clear pos; pos=find(omni(:,1)<smjd); omni(pos,:)=[];
clear pos; pos=find(omni(:,1)>fmjd); omni(pos,:)=[];

%process the spacecraft data
%============================
%tranform from GSE to RTN
[omni(:,2), omni(:,3), omni(:,4)]=GSE2heliograph(omni(:,1), omni(:,2), omni(:,3), omni(:,4));
%compute |B|
omni(:,15)=sqrt(omni(:,2).*omni(:,2) + omni(:,3).*omni(:,3) + omni(:,4).*omni(:,4));



%create a "geoeffectiveness" parameter
alpha=0.5;
theta=atan2(-omni(:,3),omni(:,4));
clear geoeff; geoeff(:,1)=omni(:,1);
geoeff(:,2)=(omni(:,8).^(2/3-alpha)).*(omni(:,15).^(2*alpha))...
    .*((-omni(:,5)).^(7/3 - 2*alpha)).*(sin(theta/2)).^4;

geoeff(:,2)=geoeff(:,2)/1e6;

%remove datagaps
clear pos; pos=find(isnan(geoeff(:,2))); geoeff(pos,:)=[]; omni(pos,:)=[];

gmax=max(geoeff(:,2)); gmin=min(geoeff(:,2));
%get the required percentile
thresh=percentile(geoeff(:,2),percent);

figure;
plot(mjd2fracyear(geoeff(:,1)),geoeff(:,2),'.');
set(gca,'YLim',[0 gmax/10]); ylabel('Geoeffectiveness');
hold on;
plot(get(gca,'XLim'),thresh*[1 1])

%create a sheath and ICME mask, add mean properties over the ICME
geoeff(:,3:4)=0;
geoeff(:,5)=NaN; %mean ICME speed
geoeff(:,6)=NaN; %mean ICME Bz
geoeff(:,7)=NaN; %mean ICME B
clear v B bz
for i=1:length(ICMEs)
    %mean ICME/sheath speed 
    clear pos;
    pos=find(geoeff(:,1)>=ICMEs(i,1) & geoeff(:,1)<=ICMEs(i,3));
    if length(nonans(abs(omni(pos,5))))>0
        v(i)=mean(nonans(abs(omni(pos,5))));
    else
        v(i)=NaN;
    end
    if length(nonans(abs(omni(pos,4))))>0
        bz(i)=mean(nonans(abs(omni(pos,4))));
    else
        bz(i)=NaN;
    end
    if length(nonans(abs(omni(pos,15))))>0
        B(i)=mean(nonans(abs(omni(pos,15))));
    else
        B(i)=NaN;
    end
    
    %sheath
    clear pos;
    pos=find(geoeff(:,1)>=ICMEs(i,1) & geoeff(:,1)<ICMEs(i,2));
    geoeff(pos,3)=1;
    
    %ICME type
    geoeff(pos,4)=ICMEs(i,4);
    geoeff(pos,5)=v(i);
    geoeff(pos,6)=bz(i);
    geoeff(pos,7)=B(i);
    
    %ICME
    clear pos;
    pos=find(geoeff(:,1)>=ICMEs(i,2) & geoeff(:,1)<=ICMEs(i,3));
    geoeff(pos,3)=2;
    
    %ICME type
    geoeff(pos,4)=ICMEs(i,4);
    geoeff(pos,5)=v(i);
    geoeff(pos,6)=bz(i);
    geoeff(pos,7)=B(i);
    
end


%find the terciles/quartiles/etc of the ICME speeds
for n=1:ncat-1
    %clear vn
    eval(['clear v',int2str(n)]);
    %vn=percentile(v,n/ncats);
    eval(['v',int2str(n),'=percentile(v,',int2str(n),'/ncat);']);
    
    %clear Bn
    eval(['clear B',int2str(n)]);
    %Bn=percentile(B,n/ncats);
    eval(['B',int2str(n),'=percentile(B,',int2str(n),'/ncat);']);
end
%%

%create cdfs over the required parameter ranges

bin_edges=[gmin:(gmax-gmin)/(nbins):gmax];
bin_centres=(bin_edges(2:end)+bin_edges(1:end-1))/2;

%find different solar wind populations
allnogap=find(~isnan(geoeff(:,2)));
noicme=find(geoeff(:,3)==0);
icmeandsheath=find(geoeff(:,3)>0);

%find the v and B ICME categories
icmeV1=find(geoeff(:,5)<=v1);
icmeB1=find(geoeff(:,7)<=B1);
for n=2:ncat-1
    nstr=int2str(n);
    nminusstr=int2str(n-1);
    %icmeBn=find(geoeff(:,7)>Bn-1 & geoeff(:,7)<=Bn);
    eval(['icmeV',nstr,'=find(geoeff(:,5)>v',nminusstr,' & geoeff(:,5)<=v',nstr',');']);
    eval(['icmeB',nstr,'=find(geoeff(:,7)>B',nminusstr,' & geoeff(:,7)<=B',nstr',');']);
end
%icmeBncat=find(geoeff(:,7)>Bn-1;
nstr=int2str(ncat);
nminusstr=int2str(ncat-1);
eval(['icmeV',nstr,'=find(geoeff(:,5)>v',nminusstr,');']);
eval(['icmeB',nstr,'=find(geoeff(:,7)>B',nminusstr,');']);

%find the VnBm categories
for i=1:ncat
    vstr=int2str(i);
    for j=1:ncat
        Bstr=int2str(j);
        eval(['icmeV',vstr,'B',Bstr,'=intersect(icmeV',vstr,', icmeB',Bstr,');']);
    end
end

clear cdf_all cdf_noicme cdf_icemandsheath
cdf_all=cumdf(geoeff(allnogap,2),bin_edges);
cdf_noicme=cumdf(geoeff(noicme,2),bin_edges);
cdf_icmeandsheath=cumdf(geoeff(icmeandsheath,2),bin_edges);


for n=1:ncat
    %cdf_icmeVn=cumdf(geoeff(icmeVn,2),bin_edges);
    nstr=int2str(n);
    eval(['cdf_icmeV',nstr,'=cumdf(geoeff(icmeV',nstr,',2),bin_edges);']);
    eval(['cdf_icmeB',nstr,'=cumdf(geoeff(icmeB',nstr,',2),bin_edges);']);
end


%find the VnBm categories
for i=1:ncat
    vstr=int2str(i);
    for j=1:ncat
        Bstr=int2str(j);
        %cdf_icmeViBj=cumdf(geoeff(icmeViBj,2),bin_edges);
        eval(['cdf_icmeV',vstr,'B',Bstr,'=cumdf(geoeff(icmeV',vstr,'B',Bstr,',2),bin_edges);'])
    end
end


%probabilities of exceeding threshold
p_all=1-percent
clear pos; pos=find(cdf_noicme(:,1)<thresh);
p_noicme=1-cdf_noicme(pos(end),2)
clear pos; pos=find(cdf_icmeandsheath(:,1)<thresh);
p_icmeandsheath=1-cdf_icmeandsheath(pos(end),2)

%probabilities for Bn and Vn categories
for n=1:ncat
    nstr=int2str(n);
    clear pos;
    %pos=find(cdf_icmeVn(:,1)<thresh);
    eval(['pos=find(cdf_icmeV',nstr,'(:,1)<thresh);']);
    %p_icmeVn=1-cdf_icmeVn(pos(end),2)
    eval(['p_icmeV',nstr,'=1-cdf_icmeV',nstr,'(pos(end),2);']);
    
    clear pos;
    eval(['pos=find(cdf_icmeB',nstr,'(:,1)<thresh);']);
    eval(['p_icmeB',nstr,'=1-cdf_icmeB',nstr,'(pos(end),2);']);
end

%probabilities for BnVm categories
for i=1:ncat
    vstr=int2str(i);
    for j=1:ncat
        Bstr=int2str(j);
        VBstr=['V',vstr,'B',Bstr];
        
        clear pos;
        %pos=find(cdf_icmeViBj(:,1)<thresh);
        eval(['pos=find(cdf_icme',VBstr,'(:,1)<thresh);']);
        %p_icmeViBj=1-cdf_icmeViBj(pos(end),2)
        eval(['p_icme',VBstr,'=1-cdf_icme',VBstr,'(pos(end),2);']);
    end
end

%%
costs=[1/ncl:1/ncl:1-1/ncl]';
cost_clim=0*ones(length(costs),1);
cost_icmes=0*ones(length(costs),1);
cost_perfect=0*ones(length(costs),1);
cost_V=0*ones(length(costs),1);
cost_B=0*ones(length(costs),1);
cost_VB=0*ones(length(costs),1);

for icost=1:length(costs)
    cost=costs(icost);
    loss=1;
    
    %number of observed intervals which exceed threshold
    nthreshexceed=length(find(geoeff(:,2)>thresh));
    %numebr of data points
    nall =length(nonans(geoeff(:,2)));
    %cost of a perfect deterministic forecast is the number of times the
    %threshold is exceeded, multiplied by the cost
    cost_perfect(icost)=cost*nthreshexceed;
    
    %cost of climatology is different if the climatological probability is
    %above or below the c/l ratio
    if (cost<=p_all)  %always take action
        %cost is the clratio at all times
        cost_clim(icost)=cost*nall;
    elseif (cost>p_all) %never take action
        %cost is simply all the missed events
        cost_clim(icost)=loss*nthreshexceed;
    end
    
    %ICME/noICME forecast
    %======================================================================
    nabovethresh_noicme=length(find(geoeff(noicme,2)>thresh));
    nbelowthresh_noicme=length(find(geoeff(noicme,2)<=thresh));
    n_noicme=nabovethresh_noicme + nbelowthresh_noicme;
    
    nabovethresh_icmeandsheath=length(find(geoeff(icmeandsheath,2)>thresh));
    nbelowthresh_icmeandsheath=length(find(geoeff(icmeandsheath,2)<=thresh));
    n_icmeandsheath=nabovethresh_icmeandsheath + nbelowthresh_icmeandsheath;
    
    %first the non-icme solar wind
    if (cost<=p_noicme) %always take action
        cost_icmes(icost)=n_noicme*cost;
    elseif(cost>p_noicme) %never take action, cost is missed events
        cost_icmes(icost)=nabovethresh_noicme*loss;
    end
    
    %then add the icme+sheath solar wind
    if (cost<=p_icmeandsheath) %always take action
        cost_icmes(icost)=cost_icmes(icost) + n_icmeandsheath*cost;
    elseif(cost>p_icmeandsheath) %never take action, cost is missed events
        cost_icmes(icost)=cost_icmes(icost) + nabovethresh_icmeandsheath*loss;
    end
       
    %ICME V and B categories 
    %======================================================================
    for n=1:ncat
        nstr=int2str(n);
        %nabovethresh_icmeV1=length(find(geoeff(icmeV1,2)>thresh));
        eval(['nabovethresh_icmeV',nstr,'=length(find(geoeff(icmeV',nstr,',2)>thresh));']);
        eval(['nabovethresh_icmeB',nstr,'=length(find(geoeff(icmeB',nstr,',2)>thresh));']);
        
        %nbelowthresh_icmeV1=length(find(geoeff(icmeV1,2)<=thresh));
        eval(['nbelowthresh_icmeV',nstr,'=length(find(geoeff(icmeV',nstr,',2)<=thresh));']);
        eval(['nbelowthresh_icmeB',nstr,'=length(find(geoeff(icmeB',nstr,',2)<=thresh));']);
        
        %n_icmeV1=nabovethresh_icmeV1 + nbelowthresh_icmeV1;
        eval(['n_icmeV',nstr,'=nabovethresh_icmeV',nstr,' + nbelowthresh_icmeV',nstr,';']);
        eval(['n_icmeB',nstr,'=nabovethresh_icmeB',nstr,' + nbelowthresh_icmeB',nstr,';']);
    end
    
    %first the non-icme solar wind
    if (cost<=p_noicme) %always take action
        cost_V(icost)=n_noicme*cost;
        cost_B(icost)=n_noicme*cost;
    elseif(cost>p_noicme) %never take action, cost is missed events
        cost_V(icost)=nabovethresh_noicme*loss;
        cost_B(icost)=nabovethresh_noicme*loss;
    end
    
    %then add the ICME categories
    for n=1:ncat
        nstr=int2str(n);
        
        eval(['if (cost<=p_icmeV',nstr,') cost_V(icost)=cost_V(icost) + n_icmeV',nstr,'*cost;'...
            'elseif(cost>p_icmeV',nstr,') cost_V(icost)=cost_V(icost) + nabovethresh_icmeV',nstr,'*loss; end']);
        eval(['if (cost<=p_icmeB',nstr,') cost_B(icost)=cost_B(icost) + n_icmeB',nstr,'*cost;'...
            'elseif(cost>p_icmeB',nstr,') cost_B(icost)=cost_B(icost) + nabovethresh_icmeB',nstr,'*loss; end']);
    end
    
    %combined ICME V and B categories
    %================================
    for i=1:ncat
        vstr=int2str(i);
        for j=1:ncat
            Bstr=int2str(j);
            VBstr=['V',vstr,'B',Bstr];
            
            %nabovethresh_icmeViBj=length(find(geoeff(icmeViBj,2)>thresh));
            eval(['nabovethresh_icme',VBstr,'=length(find(geoeff(icme',VBstr,',2)>thresh));']);
            
            %nbelowthresh_icmeViBj=length(find(geoeff(icmeViBj,2)<=thresh));
            eval(['nbelowthresh_icme',VBstr,'=length(find(geoeff(icme',VBstr,',2)<=thresh));']);
            
            %n_icmeViBj=nabovethresh_icmeViBj + nbelowthresh_icmeViBj;
            eval(['n_icme',VBstr,'=nabovethresh_icme',VBstr,' + nbelowthresh_icme',VBstr,';']);
        end
    end
    
    %first the non-icme solar wind
    if (cost<=p_noicme) %always take action
        cost_VB(icost)=n_noicme*cost;
    elseif(cost>p_noicme) %never take action, cost is missed events
        cost_VB(icost)=nabovethresh_noicme*loss;
    end
    
    %then add the ICME categories
    for i=1:ncat
        for j=1:ncat
            VBstr=['V',int2str(i),'B',int2str(j)];
            
            eval(['if (cost<=p_icme',VBstr,') cost_VB(icost)=cost_VB(icost) + n_icme',VBstr,'*cost;'...
                'elseif(cost>p_icme',VBstr,') cost_VB(icost)=cost_VB(icost) + nabovethresh_icme',VBstr,'*loss; end']);
       end
    end
end

%%
figure;
ploterrorband(costs,costs*0,100*(cost_clim-cost_icmes)./(cost_clim-cost_perfect),[0.7 0.7 0.7]);
hold on; box on;
plot(costs,100*(cost_clim-cost_V)./(cost_clim-cost_perfect),'b');
plot(costs,100*(cost_clim-cost_B)./(cost_clim-cost_perfect),'r');
plot(costs,100*(cost_clim-cost_VB)./(cost_clim-cost_perfect),'k');
set(gca,'XScale','log');
ylabel('Forecast value (%)'); xlabel('C/L: Relative cost of taking mitigating action');
legend('ICME arrival time only','ICME V','ICME B','ICME V and B');
