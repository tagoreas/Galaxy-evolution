const static double kinv5[625]={0.50938573318431102969008756056795776,-0.201999054135754018879300868242938103,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,-0.171617696998178203155779765713693344,0.043888547718481721701661965177656815,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,0.04512844678562585311935256370572501,0.0232213372836719097279322344720453274,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,0.0156197644369657136399942587003296616,-0.0226798740595273063381248993584034593,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.33082548660836409000089323769691621,-0.029076871581563580834947425304830934,-0.055921078150567034253562327924698337,-0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,-0.201999054135754018879300868242938103,0.47606736897589650250108005230538473,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,0.04120416580489373123761450639677008,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.047489769608872611307421749937661943,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,-0.0212802077523094021808359624598092963,0.0259329067334364259339882294506009867,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.12929499014187912021991128618216064,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.47606736897589650250108005230538473,-0.208705124825054220766464055078349331,0.04120416580489373123761450639677008,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.047489769608872611307421749937661943,0.0156197644369657136399942587003296616,-0.0212802077523094021808359624598092963,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.12929499014187912021991128618216064,-0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,-0.208705124825054220766464055078349331,0.46942845198569621631777061736028368,0.0086024374158155004863305786489671871,0.04512844678562585311935256370572501,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,0.048348773637094730318427987761544205,-0.209392499754775496425226410819288084,-0.0226798740595273063381248993584034593,0.0259329067334364259339882294506009867,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.0131460493244299902690367595813393764,-0.126186637712880924912638402676012118,-0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,-0.171617696998178203155779765713693344,0.04120416580489373123761450639677008,0.04120416580489373123761450639677008,0.0086024374158155004863305786489671871,0.65082689318597309796008104638122202,-0.171617696998178203155779765713693344,-0.171617696998178203155779765713693344,0.04120416580489373123761450639677008,0.04120416580489373123761450639677008,0.0086024374158155004863305786489671871,0.0086024374158155004863305786489671871,-0.031439860068336528554438658839402338,-0.171617696998178203155779765713693344,0.04120416580489373123761450639677008,0.04120416580489373123761450639677008,0.0086024374158155004863305786489671871,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,-0.053637747404387036538383169415883882,0.84697733466704012906162081167735024,0.128773269508504881924166907389650801,-0.128773269508504881924166907389650801,0.128773269508504881924166907389650801,0.128773269508504881924166907389650801,0.043888547718481721701661965177656815,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,0.04512844678562585311935256370572501,-0.171617696998178203155779765713693344,0.50938573318431102969008756056795776,0.043888547718481721701661965177656815,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,0.04512844678562585311935256370572501,-0.0226798740595273063381248993584034593,0.0232213372836719097279322344720453274,0.043888547718481721701661965177656815,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,0.04512844678562585311935256370572501,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,0.245827536876233474912383484467386944,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,0.055921078150567034253562327924698337,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,-0.171617696998178203155779765713693344,0.043888547718481721701661965177656815,0.50938573318431102969008756056795776,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,-0.0226798740595273063381248993584034593,0.04512844678562585311935256370572501,0.0232213372836719097279322344720453274,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,0.245827536876233474912383484467386944,-0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,0.055921078150567034253562327924698337,-0.029076871581563580834947425304830934,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,0.04120416580489373123761450639677008,-0.201999054135754018879300868242938103,-0.201999054135754018879300868242938103,0.47606736897589650250108005230538473,-0.0212802077523094021808359624598092963,0.0259329067334364259339882294506009867,0.0259329067334364259339882294506009867,0.047489769608872611307421749937661943,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,-0.044838547173864441761024471867525378,-0.094867545488743565676208572421196756,0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.079265991826999996304727185628489262,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,0.04120416580489373123761450639677008,0.0156197644369657136399942587003296616,0.0156197644369657136399942587003296616,-0.0212802077523094021808359624598092963,0.47606736897589650250108005230538473,-0.208705124825054220766464055078349331,-0.208705124825054220766464055078349331,0.047489769608872611307421749937661943,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.12929499014187912021991128618216064,0.079265991826999996304727185628489262,-0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,-0.094867545488743565676208572421196756,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,0.0086024374158155004863305786489671871,0.04512844678562585311935256370572501,-0.0226798740595273063381248993584034593,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.46942845198569621631777061736028368,0.048348773637094730318427987761544205,-0.209392499754775496425226410819288084,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,-0.0131460493244299902690367595813393764,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,-0.126186637712880924912638402676012118,0.090815707806348549388871629006537109,0.04512844678562585311935256370572501,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.048348773637094730318427987761544205,0.0086024374158155004863305786489671871,-0.0226798740595273063381248993584034593,0.04512844678562585311935256370572501,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.048348773637094730318427987761544205,0.46942845198569621631777061736028368,-0.209392499754775496425226410819288084,0.04512844678562585311935256370572501,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.048348773637094730318427987761544205,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,-0.0131460493244299902690367595813393764,0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,-0.126186637712880924912638402676012118,0.0232213372836719097279322344720453274,0.047489769608872611307421749937661943,0.047489769608872611307421749937661943,-0.209392499754775496425226410819288084,-0.031439860068336528554438658839402338,0.0232213372836719097279322344720453274,0.0232213372836719097279322344720453274,0.047489769608872611307421749937661943,0.047489769608872611307421749937661943,-0.209392499754775496425226410819288084,-0.209392499754775496425226410819288084,0.46822666098313243800020237145874682,0.0232213372836719097279322344720453274,0.047489769608872611307421749937661943,0.047489769608872611307421749937661943,-0.209392499754775496425226410819288084,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,0.0057398078290956923747206232859137222,-0.13307106891710070463141224853511895,0.097843568762155153181163906469981746,-0.097843568762155153181163906469981746,0.097843568762155153181163906469981746,0.097843568762155153181163906469981746,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,0.0156197644369657136399942587003296616,-0.0226798740595273063381248993584034593,-0.171617696998178203155779765713693344,0.043888547718481721701661965177656815,0.043888547718481721701661965177656815,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,0.04512844678562585311935256370572501,0.0232213372836719097279322344720453274,0.50938573318431102969008756056795776,-0.201999054135754018879300868242938103,-0.201999054135754018879300868242938103,0.04512844678562585311935256370572501,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,0.245827536876233474912383484467386944,0.055921078150567034253562327924698337,0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,-0.0212802077523094021808359624598092963,0.0259329067334364259339882294506009867,0.04120416580489373123761450639677008,-0.201999054135754018879300868242938103,0.0156197644369657136399942587003296616,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.047489769608872611307421749937661943,-0.201999054135754018879300868242938103,0.47606736897589650250108005230538473,0.047056700500673124240131459797551201,-0.208705124825054220766464055078349331,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,-0.044838547173864441761024471867525378,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,-0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.0156197644369657136399942587003296616,-0.0212802077523094021808359624598092963,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,0.04120416580489373123761450639677008,0.0156197644369657136399942587003296616,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.047056700500673124240131459797551201,0.0259329067334364259339882294506009867,-0.208705124825054220766464055078349331,0.047489769608872611307421749937661943,-0.201999054135754018879300868242938103,0.047056700500673124240131459797551201,0.47606736897589650250108005230538473,-0.208705124825054220766464055078349331,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.044838547173864441761024471867525378,0.079265991826999996304727185628489262,0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,-0.0226798740595273063381248993584034593,0.0259329067334364259339882294506009867,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,0.0086024374158155004863305786489671871,0.04512844678562585311935256370572501,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,0.0259329067334364259339882294506009867,0.048348773637094730318427987761544205,0.048348773637094730318427987761544205,-0.209392499754775496425226410819288084,0.04512844678562585311935256370572501,-0.208705124825054220766464055078349331,-0.208705124825054220766464055078349331,0.46942845198569621631777061736028368,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,-0.230148394843659464570546791263888604,0.090815707806348549388871629006537109,0.126186637712880924912638402676012118,0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,-0.0092273327027093340061026585412559696,-0.053637747404387036538383169415883882,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0057398078290956923747206232859137222,0.061995248417598442138174881428363804,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.104066129739236438129786732047189479,-0.3359902541872321754844773217002408,-0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.061995248417598442138174881428363804,0.0241941425336893278105239728891133004,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,-0.053637747404387036538383169415883882,-0.036072620747944560269137923139530078,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,0.0103846033740640180498128650637113438,0.0057398078290956923747206232859137222,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,-0.0208965840633061976516446304038626469,-0.0092273327027093340061026585412559696,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.45379875985638922477789543103573895,0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,-0.053637747404387036538383169415883882,-0.036072620747944560269137923139530078,0.061995248417598442138174881428363804,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,-0.0092273327027093340061026585412559696,0.0103846033740640180498128650637113438,0.0057398078290956923747206232859137222,-0.036072620747944560269137923139530078,-0.0208965840633061976516446304038626469,0.0241941425336893278105239728891133004,0.0103846033740640180498128650637113438,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.036389964466353233293659516186061631,-0.104066129739236438129786732047189479,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,-0.036072620747944560269137923139530078,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0103846033740640180498128650637113438,-0.053637747404387036538383169415883882,0.061995248417598442138174881428363804,-0.036072620747944560269137923139530078,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0103846033740640180498128650637113438,-0.0092273327027093340061026585412559696,0.0057398078290956923747206232859137222,-0.036072620747944560269137923139530078,0.0241941425336893278105239728891133004,-0.0208965840633061976516446304038626469,0.0103846033740640180498128650637113438,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,-0.036389964466353233293659516186061631,0.171471293969954172093906018572751458,-0.104066129739236438129786732047189479,0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,0.33082548660836409000089323769691621,0.12929499014187912021991128618216064,0.12929499014187912021991128618216064,-0.0131460493244299902690367595813393764,0.84697733466704012906162081167735024,0.245827536876233474912383484467386944,0.245827536876233474912383484467386944,-0.044838547173864441761024471867525378,0.12929499014187912021991128618216064,-0.0131460493244299902690367595813393764,-0.0131460493244299902690367595813393764,-0.13307106891710070463141224853511895,0.245827536876233474912383484467386944,-0.044838547173864441761024471867525378,-0.044838547173864441761024471867525378,-0.230148394843659464570546791263888604,-0.104066129739236438129786732047189479,-0.45379875985638922477789543103573895,-0.104066129739236438129786732047189479,-0.104066129739236438129786732047189479,4.1396342447084208887240462702559604,-0.75818538474656752840816530442692083,-2.86446276038217932737512359936226984,-0.75818538474656752840816530442692083,-0.75818538474656752840816530442692083,-0.029076871581563580834947425304830934,-0.094867545488743565676208572421196756,-0.094867545488743565676208572421196756,-0.126186637712880924912638402676012118,0.128773269508504881924166907389650801,-0.029076871581563580834947425304830934,-0.029076871581563580834947425304830934,-0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.090815707806348549388871629006537109,0.090815707806348549388871629006537109,0.097843568762155153181163906469981746,0.055921078150567034253562327924698337,0.079265991826999996304727185628489262,0.079265991826999996304727185628489262,0.090815707806348549388871629006537109,-0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.75818538474656752840816530442692083,3.824188374201590623804406888784741,-0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,-0.055921078150567034253562327924698337,-0.079265991826999996304727185628489262,-0.079265991826999996304727185628489262,-0.090815707806348549388871629006537109,-0.128773269508504881924166907389650801,0.029076871581563580834947425304830934,0.029076871581563580834947425304830934,0.094867545488743565676208572421196756,-0.079265991826999996304727185628489262,-0.090815707806348549388871629006537109,-0.090815707806348549388871629006537109,-0.097843568762155153181163906469981746,0.029076871581563580834947425304830934,0.094867545488743565676208572421196756,0.094867545488743565676208572421196756,0.126186637712880924912638402676012118,-0.0137423759299206111636313772883086637,0.3359902541872321754844773217002408,-0.0137423759299206111636313772883086637,-0.0137423759299206111636313772883086637,-2.86446276038217932737512359936226984,-0.201540229072843768021117984995550364,3.824188374201590623804406888784741,-0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,-0.029076871581563580834947425304830934,-0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.090815707806348549388871629006537109,0.128773269508504881924166907389650801,-0.029076871581563580834947425304830934,0.055921078150567034253562327924698337,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,-0.126186637712880924912638402676012118,0.090815707806348549388871629006537109,0.097843568762155153181163906469981746,-0.029076871581563580834947425304830934,-0.094867545488743565676208572421196756,0.079265991826999996304727185628489262,0.090815707806348549388871629006537109,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,0.0137423759299206111636313772883086637,-0.75818538474656752840816530442692083,0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,3.824188374201590623804406888784741,0.201540229072843768021117984995550364,-0.029076871581563580834947425304830934,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.090815707806348549388871629006537109,0.128773269508504881924166907389650801,0.055921078150567034253562327924698337,-0.029076871581563580834947425304830934,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.090815707806348549388871629006537109,-0.126186637712880924912638402676012118,0.097843568762155153181163906469981746,-0.029076871581563580834947425304830934,0.079265991826999996304727185628489262,-0.094867545488743565676208572421196756,0.090815707806348549388871629006537109,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,0.0137423759299206111636313772883086637,-0.3359902541872321754844773217002408,-0.75818538474656752840816530442692083,0.201540229072843768021117984995550364,-0.201540229072843768021117984995550364,0.201540229072843768021117984995550364,3.824188374201590623804406888784741};