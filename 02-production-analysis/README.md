## The Production Data

The production data is in the file `sp-production.csv`, with these columns and variables:

variable|type|explanation
:-:|:-:|:-:
`score`|dependent variable| whether participant produced the target structure (1) or not (0)
`priming_amount` | independent variable | the total number of trials where participant produced a stranded preposition after a prime trial during both alignment sessions
`group` | independent variable | whether participant was in `control` condition (no alignment sessions) or `exp` condition (completed alignment sessions)
`modality` | independent variable | whether participant was in an FTF or SCMC context
`time` | independent variable | test order: pretest, immediate posttest, delayed posttest
`wmc` | independent variable | particpant's working memory capacity score 
`rec_pre` | independent variable | participant's receptive knowledge score (GJT)
`cloze` | independent variable | participant's proficiency score (cloze test)
`test` | control variable | test version (1, 2, or 3, for counterbalancing)
`trial_order` |control variable| order of the questions within any one production test session
`subject` | random effect | random intercept fit for each subject
`verb`| random effect | random intercept fit for each test question (each question had a unique verb)
