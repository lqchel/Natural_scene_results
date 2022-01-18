function [out] = AUC_lme(data, hypo)

if hypo == 1 % do not take congruence into account
    data_table = table(categorical(data(data(:,2)~= -1,2)),data(data(:,2)~= -1,3),...
     data(data(:,2)~= -1,1), 'VariableNames',{'eccentricity','AUC','subjects'});
    lm_full = fitlme(data_table,'AUC~eccentricity + (eccentricity|subjects)+(1|subjects)');
    lm_reduced = fitlme(data_table, 'AUC ~ 1 + (eccentricity|subjects) + (1|subjects)');
    out = compare(lm_reduced,lm_full);
    
elseif hypo == 2

        data_table = table(categorical(data(data(:,3)~= -1,3)),data(data(:,3)~= -1,4),...
            data(data(:,3)~= -1,1),categorical(data(data(:,3)~= -1,2)),'VariableNames',{'eccentricity','AUC','subjects','condition'});
        lme_full = fitlme(data_table,'AUC~eccentricity*condition + (eccentricity|subjects) + (condition|subjects) + (1|subjects)');
        lme_ecc = fitlme(data_table,'AUC~eccentricity:condition + condition + (eccentricity|subjects) + (condition|subjects)+ (1|subjects)');
        lme_interaction = fitlme(data_table,'AUC~eccentricity + condition + (eccentricity|subjects) + (condition|subjects)+ (1|subjects)');
        lme_condition = fitlme(data_table,'AUC~eccentricity + eccentricity:condition + (eccentricity|subjects) + (condition|subjects)+ (1|subjects)');
        
        interaction = compare(lme_interaction, lme_full)
        ecc = compare(lme_ecc,lme_full)
        condition = compare(lme_condition,lme_full)

        out = {ecc, interaction, condition}
    
    
end
end