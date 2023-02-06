library(tidyverse)
library(data.table)
library(lubridate)
library(lme4)
library(broom.mixed)

vega_cohort <- read.csv('~/vega-ars/3_output/VEGA_COHORT.csv')

create_forest_plot <- function (df, title = NA, subtitle = NA)
{
  df$lower <- df$estimate - 1.96 * df$std.error
  df$upper <- df$estimate + 1.96 * df$std.error
  df$exp_coef <- exp(df$estimate)
  df$exp_std.error <- exp(df$std.error)
  df$exp_lower <- exp(df$lower)
  df$exp_upper <- exp(df$upper)
  df$term <- factor(df$term, levels = rev(df$term))
  df <- dplyr::arrange(df, estimate)
  df$term <- factor(df$term, levels = rev(df$term))
  forest_plot_exp <- ggplot2::ggplot(data = df,
                                     ggplot2::aes(x = term,
                                                  y = exp_coef,
                                                  ymin = exp_lower,
                                                  ymax = exp_upper)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept = 1, lty = 2) +
    ggplot2::coord_flip() +
    ggplot2::ylab("Odds Ratio (95% CI)") +
    ggplot2::xlab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 15))

  df$exp_coef <- round(df$exp_coef, 2)
  df$exp_lower <- round(df$exp_lower, 2)
  df$exp_upper <- round(df$exp_upper, 2)
  df$p.value <- round(df$p.value, 3)
  if (!is.na(title)) {
    forest_plot_exp <- forest_plot_exp + ggplot2::ggtitle(title)
  }
  if (!is.na(subtitle)) {
    forest_plot_exp <- forest_plot_exp + ggplot2::labs(title = title,
                                                       subtitle = subtitle)
  }
  return(forest_plot_exp)
}

### BEGIN FIGURE 3 ###
analyze_secondary_outcome <- function(df, secondary_outcome, name) {
  set.seed(1234); summary(secondary_outcome_model <- glmer(
    as.formula(paste0(secondary_outcome, '~
                      ASSIGNED_NOREPI +
                      PAT_AGE_AT_ENC +
                      GENDER +
                      ETHNICITY +
                      SURGICAL_SERVICE +
                      CASE_DURATION +
                      RACE_ASIAN +
                      RACE_BLACK +
                      RACE_WHITE +
                      RACE_NATIVE +
                      RACE_PACIFIC +
                      RACE_OTHER +
                      PREOPERATIVE_HYPERTENSION_FLAG +
                      PREOPERATIVE_ACE_ARB_FLAG +
                      PREOPERATIVE_DIABETES_FLAG +
                      PREOPERATIVE_CHF_FLAG +
                      PREOPERATIVE_CAD_FLAG +
                      HYPOTENSION_DURATION_MINUTES +
                      ne_bolus_count_zeroed +
                      pe_bolus_count_zeroed +
                      ephedrine_bolus_count_zeroed +
                      multiple_pressor_infusions +
                      (1 | CENTER_NAME)')),
    data=df,
    na.action=na.exclude,
    nAGQ = 0,
    family=binomial(link='logit'),
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  ))

  secondary_outcome_broom <- broom.mixed::tidy(secondary_outcome_model) %>%
    filter(term == 'ASSIGNED_NOREPITRUE') %>%
    mutate(term = name)

  secondary_outcome_df <<- rbind(secondary_outcome_df, secondary_outcome_broom)
  secondary_outcome_df
}

secondary_outcome_df <- data.frame()
analyze_secondary_outcome(vega_cohort, 'ANY_AKI',  'AKI')
analyze_secondary_outcome(vega_cohort, 'MORTALITY',  '30-Day Mortality')
analyze_secondary_outcome(vega_cohort, 'MINS',  'MINS')
analyze_secondary_outcome(vega_cohort, 'ADVERSE_CARDIO_RENAL_EVENTS',  'Adverse Cardiorenal Events')
analyze_secondary_outcome(vega_cohort, 'readmission',  'Hospital Readmission via ED')
analyze_secondary_outcome(vega_cohort, 'SEVERE_AKI_BOOL',  'Severe AKI')
secondary_outcome_df

figure_3 <- create_forest_plot(secondary_outcome_df,
                               title='Odds Ratios of Secondary Outcomes',
                               subtitle = 'Patients Randomized to NE compared to PE')
### END FIGURE 3 ###

### BEGIN eFIGURE 2 ###
set.seed(1024)
aki_model <- glmer(
  ANY_AKI ~
    ASSIGNED_NOREPI +
    PAT_AGE_AT_ENC +
    GENDER +
    ETHNICITY +
    SURGICAL_SERVICE +
    CASE_DURATION +
    RACE_ASIAN +
    RACE_BLACK +
    RACE_WHITE +
    RACE_NATIVE +
    RACE_PACIFIC +
    RACE_OTHER +
    PREOPERATIVE_HYPERTENSION_FLAG +
    PREOPERATIVE_ACE_ARB_FLAG +
    PREOPERATIVE_DIABETES_FLAG +
    PREOPERATIVE_CHF_FLAG +
    PREOPERATIVE_CAD_FLAG +
    HYPOTENSION_DURATION_MINUTES +
    pe_bolus_count_zeroed +
    ne_bolus_count_zeroed +
    ephedrine_bolus_count_zeroed +
    multiple_pressor_infusions +
    (1 | CENTER_NAME),
  data = vega_cohort,
  nAGQ = 0,
  family=binomial(link='logit'),
  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
)

aki_model_df <- broom.mixed::tidy(aki_model) %>%
  filter((p.value <= 0.05 | term == 'ASSIGNED_NOREPITRUE'), term != '(Intercept)') %>%
  mutate(
    term = ifelse(term == 'ASSIGNED_NOREPITRUE', 'Assigned to NE', term),
    term = ifelse(term == 'SURGICAL_SERVICEOrtho', 'Orthopedic Surgery', term),
    term = ifelse(term == 'PREOPERATIVE_HYPERTENSION_FLAGY', 'Preoperative Hypertension', term),
    term = ifelse(term == 'PREOPERATIVE_CHF_FLAGY', 'Heart Failure', term),
    term = ifelse(term == 'SURGICAL_SERVICEOther', '"Other" Surgical Service', term),
    term = ifelse(term == 'RACE_PACIFICTrue', 'Native Hawiian/Pacific Islander', term),
    term = ifelse(term == 'CASE_DURATION', 'Case Duration', term),

  )

efigure_2 <- create_forest_plot(aki_model_df, title='Odds Ratios of AKI')
### END eFIGURE 2 ###

### BEGIN eFIGURE 3 ###
sub_age <- vega_cohort %>% filter(PAT_AGE_AT_ENC >= 65)
sub_htn <- vega_cohort %>% filter(PREOPERATIVE_HYPERTENSION_FLAG == 'Y')
sub_dm <- vega_cohort %>% filter(PREOPERATIVE_DIABETES_FLAG== 'Y')
sub_chf <- vega_cohort %>% filter(PREOPERATIVE_CHF_FLAG == 'Y')
sub_acearb <- vega_cohort %>% filter(PREOPERATIVE_ACE_ARB_FLAG == 'Y')
sub_beta <- vega_cohort %>% filter(PREOPERATIVE_BETA_BLOCKER_FLAG== 'Y')
sub_neuro <- vega_cohort %>% filter(SURGICAL_SERVICE== 'Neuro')
sub_abdpelvic <- vega_cohort %>% filter(SURGICAL_SERVICE == 'Abdominal/Pelvic')
sub_ebl <- vega_cohort %>% filter(EBL_ML >= 300)

nsub_age <- vega_cohort %>% filter(PAT_AGE_AT_ENC < 65)
nsub_htn <- vega_cohort %>% filter(PREOPERATIVE_HYPERTENSION_FLAG == 'N')
nsub_dm <- vega_cohort %>% filter(PREOPERATIVE_DIABETES_FLAG== 'N')
nsub_chf <- vega_cohort %>% filter(PREOPERATIVE_CHF_FLAG == 'N')
nsub_acearb <- vega_cohort %>% filter(PREOPERATIVE_ACE_ARB_FLAG == 'N')
nsub_beta <- vega_cohort %>% filter(PREOPERATIVE_BETA_BLOCKER_FLAG== 'N')
nsub_neuro <- vega_cohort %>% filter(SURGICAL_SERVICE != 'Neuro')
nsub_abdpelvic <- vega_cohort %>% filter(SURGICAL_SERVICE != 'Abdominal/Pelvic')
nsub_ebl <- vega_cohort %>% filter(EBL_ML < 300)

analyze_subgroup <- function(df, outcome, name) {
  set.seed(1234); summary(sub_model <- glmer(
    as.formula(paste0(outcome, '~ ASSIGNED_NOREPI + (1 | CENTER_NAME)')),
    data=df,
    na.action=na.exclude,
    family=binomial(link='logit'),
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  ))

  sub_broom <- broom.mixed::tidy(sub_model) %>%
    filter(term == 'ASSIGNED_NOREPITRUE') %>%
    mutate(term = name)

  sub_df <<- rbind(sub_df, sub_broom)
  sub_df
}

# AKI
# Not in subgroup analysis
sub_df <- data.frame()
analyze_subgroup(nsub_age, 'ANY_AKI',  'Age ≥ 65')
analyze_subgroup(nsub_htn, 'ANY_AKI',  'Hypertension')
analyze_subgroup(nsub_dm, 'ANY_AKI',  'Diabetes')
analyze_subgroup(nsub_chf, 'ANY_AKI',  'Heart Failure')
analyze_subgroup(nsub_acearb, 'ANY_AKI',  'ACEi and/or ARB')
analyze_subgroup(nsub_beta, 'ANY_AKI',  'Beta Blockers')
analyze_subgroup(nsub_neuro, 'ANY_AKI',  'Neurosurgery')
analyze_subgroup(nsub_abdpelvic, 'ANY_AKI',  'Abdominal/Pelvic Surgery')
analyze_subgroup(nsub_ebl, 'ANY_AKI',  'EBL ≥ 300 mL')
nsub_df <- sub_df %>%
  mutate(group='Not in Subgroup')

# In subgroup analysis
sub_df <- data.frame()
analyze_subgroup(sub_age, 'ANY_AKI',  'Age ≥ 65')
analyze_subgroup(sub_htn, 'ANY_AKI',  'Hypertension')
analyze_subgroup(sub_dm, 'ANY_AKI',  'Diabetes')
analyze_subgroup(sub_chf, 'ANY_AKI',  'Heart Failure')
analyze_subgroup(sub_acearb, 'ANY_AKI',  'ACEi and/or ARB')
analyze_subgroup(sub_beta, 'ANY_AKI',  'Beta Blockers')
analyze_subgroup(sub_neuro, 'ANY_AKI',  'Neurosurgery')
analyze_subgroup(sub_abdpelvic, 'ANY_AKI',  'Abdominal/Pelvic Surgery')
analyze_subgroup(sub_ebl, 'ANY_AKI',  'EBL ≥ 300 mL')
aki_df <- sub_df %>%
  mutate(group='In Subgroup') %>%
  rbind(nsub_df)

aki_df$lower <- aki_df$estimate - 1.96 * aki_df$std.error
aki_df$upper <- aki_df$estimate + 1.96 * aki_df$std.error
aki_df$exp_coef <- exp(aki_df$estimate)
aki_df$exp_std.error <- exp(aki_df$std.error)
aki_df$exp_lower <- exp(aki_df$lower)
aki_df$exp_upper <- exp(aki_df$upper)
# Arrange the df in order of OR, use this order to order the factors in the main df
sub_df <- dplyr::arrange(sub_df, estimate)
aki_df$term <- factor(aki_df$term, levels = rev(sub_df$term))
forest_plot_colors <- c('#000000', '#CCCCCC')
efigure_3 <- ggplot2::ggplot(
  data = aki_df,
  aes(x = factor(term, level=rev(sub_df$term)), y = exp_coef, ymin = exp_lower, ymax = exp_upper, col=group, fill=group)) +
  geom_pointrange(position=position_dodge(width=-0.75)) +
  scale_fill_manual(values=forest_plot_colors) +
  scale_color_manual(values=forest_plot_colors) +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  ylab("Odds Ratio (95% CI)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.y = ggplot2::element_text(size = 15))  +
  labs(title='Odds Ratio of AKI within Subgroups',
       subtitle='Patients randomized to NE compared to PE',
  )
### END eFIGURE 3 ###
