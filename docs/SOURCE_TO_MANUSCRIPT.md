# Scientific output map

The following tables and figures are computed from the supplied inputs. Narrative reporting materials and upstream metadata are maintained separately.

| Object | Producer | Computed output or boundary |
|---|---|---|
| Main Table 1 | R/modules/descriptive.R;validation.R | baseline_descriptive_by_cohort.csv;baseline_group_tests_continuous.csv;baseline_group_tests_categorical_seeded.csv |
| Main Table 2 | R/modules/validation.R | primary_validation.csv |
| Main Figure 1 | R/render_figures.R | primary_1d_window_details_figure_ready.csv -> primary_1d_landscape_six_panel.pdf |
| Main Figure 2 | R/render_figures.R | primary_2d_cell_details.csv -> sips_n_brainage_local_auc_surface.pdf |
| Main Figure 4 | R/render_figures.R | canonical_linear_moderation_curve.csv -> canonical_linear_sips_n_moderation.pdf |
| Supplement Figure 1 | R/render_figures.R | external_auc_random_effects_forest_data.csv -> external_auc_random_effects_forest.pdf |
| Supplement Figure 2 | R/render_figures.R | external_original_vs_apparent_calibration_curve.csv -> external_calibration_original_vs_apparent.pdf |
| Supplement Figure 4 | R/render_figures.R | component_selected_statistic_key.csv -> component_selected_statistic_heatmap.pdf |
| Supplement Figure 8 | R/render_figures.R | gam_sips_n_uncertainty_segmented.csv -> gam_sips_n_uncertainty_supported_regions.pdf |
| Main Figure 3 | python/yeo_overlap.py;python/figure3.py | original signed MRIcroGL source_panels;yeo_recorded_full_overlap.csv;yeo_radar_display.csv |
| Supplement Figure 3 | python/aal3_overlap.py | original AAL3 .mat/.xlsx -> aal3_recorded_full_overlap.csv;aal3_display.csv;supplementary_figure3_aal3.pdf |
| Supplement Table 1 | R/modules/calibration.R;cohort_calibration.R;validation.R | external_calibration_original_vs_apparent_bootstrap.csv;calibration_audit.csv;source_calibration.csv;iecv_training_prevalence_constant.csv;sex_stratified_performance.csv |
| Supplement Table 2 | R/modules/correlations.R | correlation_family_S2_original_21.csv |
| Supplement Table 3 | R/modules/correlations.R | correlation_family_S3_mri_risk_score_7.csv |
| Supplement Table 4 | R/modules/components.R | score_component_correlations.csv;subject_level_component_moderator_correlations.csv |
| Supplement Table 5 | R/modules/components.R | table_S5_component_landscape_associations.csv;component_selected_statistic_key.csv |
| Supplement Table 6 | R/modules/moderation.R | canonical_glm_coefficients.csv;glm_hierarchy.csv;glm_hierarchy_diagnostics.csv;glm_simple_slopes.csv |
| Supplement Table 7 | R/modules/moderation.R;gam.R | canonical_glm_coefficients.csv;glm_hierarchy.csv;gam_bounded_diagnostic_sensitivity.csv;gam_fit_diagnostics.csv;gam_concurvity_canonical.csv;gam_k_diagnostics_canonical.csv;gam_k_diagnostic_row_order_sensitivity.csv;gam_bernoulli_simulation_diagnostics.csv |
| Supplement Table 12 | R/modules/landscapes.R;outputs.R | global_1d_landscape_tests.csv;global_2d_surface_tests.csv;primary_1d_window_details.csv;primary_2d_cell_details.csv;analysis_populations.csv |
| Supplement Table 14 | R/modules/descriptive.R | brainage_external_performance.csv;brainage_identity_checks.csv;baseline_descriptive_by_cohort.csv |
| Supplement Table 15 | R/modules/descriptive.R;validation.R | missingness_by_cohort.csv;baseline_descriptive_by_cohort.csv;baseline_group_tests_continuous.csv;baseline_group_tests_categorical_seeded.csv |
| Supplement Table 16 | R/modules/landscapes.R | joint_width_step_sensitivity_1d.csv;joint_width_step_sensitivity_selected_surface.csv;landscape_population_sensitivity_1d_all_six.csv;window_support_sensitivity_1d_all_six.csv;highlighted_surface_population_support_sensitivity.csv |
| Supplement Table 17 | R/modules/loco.R | loco_incremental_value.csv;loco_paired_comparison.csv;loco_training_scaling.csv;loco_scaling_invariance.csv |
| Supplement Table 18 | R/modules/outputs.R | correction_registry.csv;uncertainty_methods_registry.csv;resampling_validity_registry.csv |
