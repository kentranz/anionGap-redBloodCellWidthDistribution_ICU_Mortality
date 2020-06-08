-- v0.2:
-- Nothing is changed, just re-format for easy reading

SELECT 
	d.subject_id
	, d.hadm_id
	, d.icustay_id
	, d.gender
	, d.age
	, d.ethnicity
	, d.admission_type
	, c."AG1"
	, c."RDW1"
	, c.first_careunit  -- General Predictors

	, s.sapsii
	, s.sapsii_prob
	, s.age_score
	, s.hr_score
	, s.sysbp_score
	, s.temp_score
	, s.pao2fio2_score
	, s.uo_score
	, s.bun_score   -- SAPS-II Scores
    
    , s.wbc_score
    , s.potassium_score
    , s.sodium_score
    , s.bicarbonate_score
    , s.bilirubin_score
    , s.gcs_score
    , s.comorbidity_score
    , s.admissiontype_score

    , c.congestive_heart_failure
    , c.cardiac_arrhythmias
    , c.valvular_disease
    , c.pulmonary_circulation
    , c.peripheral_vascular 	    -- Elixhauser Scores
    
    , c.hypertension
    , c.paralysis
    , c.other_neurological
    , c.chronic_pulmonary
    , c.diabetes_uncomplicated
    , c.diabetes_complicated
    , c.hypothyroidism
    , c.renal_failure
    , c.liver_disease
    , c.peptic_ulcer
    , c.aids
    , c.lymphoma
    , c.metastatic_cancer
    , c.solid_tumor
    , c.rheumatoid_arthritis
    , c.coagulopathy
    , c.obesity
    , c.weight_loss
    , c.fluid_electrolyte
    , c.blood_loss_anemia
    , c.deficiency_anemias
    , c.alcohol_abuse
    , c.drug_abuse
    , c.psychoses
    , c.depression

    , d.los_icu
    , d.los_hospital
    , d.hospital_expire_flag --Outcomes
	, (CASE WHEN p.dod < c.outtime THEN 1 ELSE 0 END) as ICU_expire_flag
	, (CASE WHEN p.dod < c.intime + interval '30' day THEN 1 ELSE 0 END) as thirty_day_mort					   				    
    , (CASE WHEN p.dod < c.intime + interval '1' year THEN 1 ELSE 0 END) AS one_year_mortality
    , ROUND( (CAST(EXTRACT(epoch FROM p.dod - c.intime)/(60*60*24) AS numeric)), 4) AS survival

	, icu.DBSOURCE

FROM mimiciii_dev."ICUSTAY_DETAIL" d
INNER JOIN mimiciii_dev."Clinical_Data" c
ON d.icustay_id = c.icustay_id

INNER JOIN mimiciii_dev.sapsii s
ON d.icustay_id = s.icustay_id

INNER JOIN mimiciii.patients p
ON d.subject_id = p.subject_id

INNER JOIN mimiciii.icustays icu
ON icu.icustay_id = c.icustay_id

where d.age>15 
and c."AG1" is not null 
and c."RDW1" is not null 
and s.sapsii is not null 
and d.icustay_seq =1 
and d.first_icu_stay='Y';
