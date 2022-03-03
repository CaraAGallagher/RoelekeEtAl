; BatSensoryNetworksModel version 0.5

; Understanding the potential importance of bat sensory networks for locating food patches in
; varying environments. Refer to the scientific publications for detailed
; documentation:

; Roeleke, M., Schlägel, U. E., Gallagher, C. A., Pufelski, J., Blohm, T., Nathan, R., Toledo, S., Jeltsch, F., & C. C. Voigt.
; Insectivorous bats form mobile sensory networks to optimize prey localization. (Submitted.)

; Model Developed by:
; Cara A. Gallagher
; Postdoctoral researcher
; University of Potsdam
; Developed as part of the BioMove project (biomove.org)

; Development started: Feb 16th, 2021
; The model was developed and tested using NetLogo version 6.0.4. Development ended August 16th, 2021.

; debug levels:
;   0   no debugging
;   1   profiler
;   2   debugging network size estimator
;   3   debugging alignment
;   4   debugging random walk

extensions [profiler]

globals [
  foraging-area                ; position of the center of the foraging area
  turn-ang-alpha-search        ; turning angle gamma distribution alpha parameter for bats which are searching
  turn-ang-lambda-search       ; turning angle gamma distribution lambda parameter for bats which are searching
  turn-ang-alpha-hunt          ; turning angle gamma distribution alpha parameter for bats which are hunting
  turn-ang-lambda-hunt         ; turning angle gamma distribution lambda parameter for bats which are hunting
  step-len-alpha-search        ; step length gamma distribution alpha parameter for bats which are searching
  step-len-lambda-search       ; step length gamma distribution lambda parameter for bats which are searching
  step-len-alpha-hunt          ; step length gamma distribution alpha parameter for bats which are hunting
  step-len-lambda-hunt         ; step length gamma distribution lambda parameter for bats which are hunting
  prey-detection-range         ; radius where prey can be detected [m]

  ; for scenarios
  run-done                     ; boolean which triggers when run is complete
  most-done                    ; boolean which triggers when 95% of bats have found food
  network-unassigned-bats      ; list to count network size which keeps track of unassigned bats
  networks                     ; list of list of all networks
  food-found-ticks-list        ; keeps track of when bats find food
  perc-search-w-consp          ; percent of searching bats with a conspecific (collected once every 8 ticks)
  perc-cells-found             ; percent of cells occupied by hunting bats (collected once every 8 ticks)
  time-to-all-food             ; how long it takes for all bats in the colony to find food
  time-to-95-food              ; how long it takes for 95% of bats to find food
  avg-netwrk                   ; average size of networks [#]
  average-network-size         ; list of average size of networks [#]
  network-sizes                ; list of network sizes for each non-foraging bat [#]
  nearest-neighbor-distance    ; list of distances to the nearest neighbor for each bat [m]
  track-list                   ; list of bats which are being tracked
  track-outputs                ; position outputs from tracked bats
  currseed                     ; seed used for simulation run
]

breed [ bats bat ]
breed [ roosts roost ]

patches-own [
  food                    ; whether or not food is found there
  patch-id                ; the number/id of the patch
  found                   ; whether this cell has been found by a bat
]

bats-own [
  step-length             ; speed [cells per tick]
  turning-angle           ; turning angle [deg]
  leave-roost-tick        ; time which bats leave the roost relative to model start [timesteps]
  conspecific             ; nearest bat to focal bat [cells]
  dist-consp              ; distance to nearest conspecific [cells]
  attract-range           ; range which attraction behavior begins [cells]
  align-range             ; range which alignment behavior is the strongest [cells]
  avoid-range             ; range which avoidance behavior is the strongest [cells]
  attract-vect            ; attraction vector [as a list of x and y values]
  align-vect              ; alignment vector [as above]
  avoid-vect              ; avoidance vector [as above]
  rw-vect                 ; correlated random walk (RW) vector biased towards foraging area [as above]
  attract-strength        ; relative weighing of attraction vector [0-1]
  align-strength          ; relative weighing of alignment vector [0-1]
  avoid-strength          ; relative weighing of avoidance vector [0-1]
  rw-strength             ; relative weighing of RW vector [0-1]
  food-found              ; whether or not a bat has found food [boolean]
  food-cell               ; cell where food was found
  hunting-consp-ticker    ; ticker to keep track of how long bats have been flying with a conspecific that is hunting [timesteps]
  move-away-ticker        ; ticker to keep track of how long bats use only RW behavior after leaving a hunting conspecific [timesteps]
  network-ID              ; index of bat's network within network list
  turn-angle-real         ; actual turning angle of bats after all movement processes have occurred [deg]
  targ-cell               ; food cell that bat targets when it has found an occupied food cell in a patch with more food around
  flying-towards-food     ; whether bat is currently flying towards a food cell

  ; for calibration
  prev-dist-c             ; previous distance to the conspecific [cells]
  diff-dist               ; difference in distances from te previous recording to the current [m]
  prev-consp              ; ID of previous conspecific
  prev-consp-hunt         ; Whether or not previous conspecific was hunting when recording occurred
]

to setup
  ca
  reset-ticks
  ; setup world dimensions
  resize-world -40 40 -40 40
  set-patch-size 6

  ; set seed for BS experiments to match maps
  ; random-seed (behaviorspace-run-number - 1)

  setup-maps

  ;; initialize globals
  ; based on empirical movements fit with gamma distributions
  set turn-ang-alpha-search 0.66
  set turn-ang-lambda-search 2.62
  set turn-ang-alpha-hunt 0.78
  set turn-ang-lambda-hunt 0.98
  set step-len-alpha-search 3.4
  set step-len-lambda-search 5.6
  set step-len-alpha-hunt 1.5
  set step-len-lambda-hunt 6.7
  set prey-detection-range 15 / 75
  set run-done false
  set most-done false
  set food-found-ticks-list []
  set perc-search-w-consp []
  set perc-cells-found []
  set average-network-size []
  set network-sizes []
  set nearest-neighbor-distance []
  if track-bats? = true [
   set track-list []
   set track-outputs []
  ]

  ; for estimating network size:
  set network-unassigned-bats []
  set networks []

  create-bats n-bats   [             ; bats spawn in the roost area
    ; select position of the roost - based on region where bats were tagged
    setxy (0 + random-float 0.1 - random-float 0.1) ((8.67 - 40) + random-float 0.1 - random-float 0.1)
    set leave-roost-tick random 37   ; used to stagger bats leaving the roost (for 8 sec timesteps this is ~5 minutes)
    set conspecific nobody           ; initialize conspecific agentset as empty
    set food-found false             ; initialize animals as not having found food
    set flying-towards-food false    ; and not flying towards and food cells
    set targ-cell 0
    set attract-vect list 0 0        ; set vectors to empty lists with two elements
    set align-vect list 0 0
    set avoid-vect list 0 0
  ]

  if track-bats? = true [ ask n-of n-track-bats bats [ set track-list lput who track-list ] ] ; assign bats to track

  ; ask bats to trace their movements
  ifelse track-bats? = true
  [ ask bats with [member? who track-list = true] [pd] ]
  [ ask n-of 3 bats [pd] ]
end

to go
  if ( debug = 1 ) and ( ticks > 3 ) [
    profiler:start
  ]

  ; end run when all bats have found food
  if not any? bats with [food-found = false] [
    set run-done true
    collect-outputs
    stop
  ]

  ask bats with [leave-roost-tick < ticks and food-found = false] [    ; bats which have left the roost go through these steps each tick:
    set-conspecific                                                    ; try to find a conspecific and determine distance to conspecific
    check-consp-hunting                                                ; check whether their conspecific is hunting and how long they have been with a hunting conspecific
  ]

  ; collect outputs for calibration step
  if dev-steps = "calibration" [
    if track-bats? = true and ( remainder ticks 4 ) = 0 [ movement-outputs ]
  ]
  ; collect outputs for evaluation step
  if dev-steps = "evaluation" [
    if track-bats? = true [ movement-outputs ]
  ]

 ; set movement speed
 ask bats [ ifelse food-found = true
   [ set step-length random-gamma step-len-alpha-hunt step-len-lambda-hunt ]
   [ set step-length random-gamma step-len-alpha-search step-len-lambda-search ]
  ]

  ask bats with [leave-roost-tick < ticks and food-found = false] [
    ; bats which have left the roost go through these steps each tick:

    ifelse conspecific != nobody [              ; if they found a conspecific, then:
      set-ranges                                ; set the various ranges used for estimating vectors based on the foraging state of the conspecific
      if null-model? = false                    ; if not running null model,
      [ attraction                              ; calculate attraction vector
        alignment                               ; calculate alignment vector
        avoidance ]                             ; then calculate avoidance vector
    ]
    [                                           ; if no conspecific was found, set non-RW vectors to zero
      set attract-vect list 0 0
      set align-vect list 0 0
      set avoid-vect list 0 0
    ]
    random-walk                                 ; calculate random walk vector
  ]

  ask bats with [leave-roost-tick < ticks and food-found = false] [    ; ask bats to actually move in a separate step (so all calculate before moving)
    bats-move                                                          ; determine resulting movement direction and move forward
  ]

  ask bats with [food-found = true] [ hunting-fly ]           ; hunting bats fly around in their food cell

  if ( ticks >= 8 ) and ( remainder ticks  8 ) = 0 [          ; collect these outputs once per every 8 ticks
    if null-model? = false [ calculate-network-size ]         ; network size estimation (and coloration)
    update-monitors                                           ; plot outputs
    collect-outputs                                           ; collect additional outputs
  ]

  tick

  ; profiling
  if ( debug = 1 ) and ( ticks > 3 ) [
    profiler:stop
    print profiler:report
    profiler:reset
  ]

end

; basic map setup
to setup-maps

  ask patches [                                              ; initialize cells with no food and a patch-id of 0
    set food false
    set patch-id 0
    set found false
  ]

  create-roosts 1 [                                          ; initialize the roost
    setxy 0 (8.67 - 40)                                      ; select position of the roost - based on batworld from Manu
    set shape "circle 2"                                     ; set shape
    set color grey
    set size 2
  ]

  if dev-steps = "calibration" [
    set no-patch 1 + random no-food-cells                    ; used to randomize the number of patches generated (in calibration)
  ]

  ; patch generation steps:
  set foraging-area patch 0 0                                ; select position of center of foraging area
  ask foraging-area [
    let patch-no no-patch                                    ; create a number of patches equal to the slider value
    while [patch-no > 0]                                     ; go through a loop to select the center of each food patch
    [
      let rad (foraging-radius - 6)
      ask one-of patches in-radius rad with [any? patches in-radius 2 with [food = true] = false] [
        set food true
        set patch-id patch-no
        set pcolor one-of base-colors
      ]
      set patch-no patch-no - 1
    ]

    let cells-no (no-food-cells - no-patch)
    let curr-patch no-patch                                  ; distribute the food cells in these different patches equal to the slider value
    while [cells-no > 0]
    [
      let ids (list curr-patch 0)
      let rad (foraging-radius - 6)
      ask one-of patches with [ patch-id = curr-patch and (any? neighbors with [food = false] and not any? neighbors with [ member? patch-id ids = false ] and any? neighbors with [distance patch 0 0 < rad]) ] [
        if any? neighbors with [ food = false and not any? neighbors with [member? patch-id ids = false ] and any? neighbors with [distance patch 0 0 < rad]] [
          ask one-of neighbors with [ food = false and not any? neighbors with [member? patch-id ids = false ] and any? neighbors with [distance patch 0 0 < rad]]  [
            set food true
            set patch-id [patch-id] of myself
            set pcolor [pcolor] of myself
            set cells-no cells-no - 1
          ]
      ]]
      set patch-no patch-no - 1
      ifelse curr-patch = 1 [ set curr-patch no-patch ][ set curr-patch curr-patch - 1 ]
    ]
  ]

end

; find and select conspecific
to set-conspecific

  let convert 75                              ; conversion for 75 m2 cell size
  let pos-conspecifics nobody                 ; initialize possible conspecifics as nobody

  ; set possible conspecifics as bats in the attraction range
  if any? other bats in-radius (attraction-range / convert) with [leave-roost-tick < ticks] [
    set pos-conspecifics other bats in-radius (attraction-range / convert) with [leave-roost-tick < ticks]
  ]

  ifelse pos-conspecifics = nobody
  [ set conspecific nobody ]                                            ; if no conspecifics found then set conspecifics to nobody
  [
    ifelse move-away-ticker > 0 [                                       ; if bats are moving away from foraging bats (from check-consp-foraging procedure)
      ifelse any? pos-conspecifics with [food-found = false]            ; check for potential conspecifics which are not currently foraging
      [
        let search-cons pos-conspecifics with [food-found = false]
        set conspecific one-of search-cons with-min [distance myself]   ; then select the closest searching bat to keep track of
        set dist-consp distance conspecific                             ; determine the distance to the conspecific to use later
      ]
      [  set conspecific nobody ]                                       ; if no conspecifics found then set conspecifics to nobody]
    ]
    [ set conspecific one-of pos-conspecifics with-min [distance myself]; if not moving away from hunting bats, then select whichever bat is closest
      set dist-consp distance conspecific  ]                            ; determine the distance to the conspecific to use later
  ]

  ; if no conspecific or if conspecific is searching reset hunting-consp-ticker (from check-consp-hunting procedure)
  if conspecific = nobody or [food-found] of conspecific = false [ set hunting-consp-ticker 0 ]

end

; check if conspecific is hunting - to prevent bats from getting stuck flying next to a conspecific that has already found food
to check-consp-hunting

  ; if currently moving away from a hunting conspecific
  if move-away-ticker > 0 [
    set move-away-ticker move-away-ticker + 1                       ; increase the timer used to keep track of how long to move away from a hunting conspecific
  ]

  if conspecific != nobody and [food-found] of conspecific = true [ ; if conspecific is hunting then use local enhancement
    if local-enhancement? = true
    [
      local-enhancement
    ]

    set hunting-consp-ticker hunting-consp-ticker + 1               ; increase the timer used to keep track of how long they've been flying next to the hunting conspecific
    if hunting-consp-ticker >= 8 * 5 [ set move-away-ticker 1 ]     ; when reaching max value, start moving away
  ]

  if move-away-ticker >= 8 * 5 [                                    ; if the timer to move away reaches its max value, reset both values
    set hunting-consp-ticker 0
    set move-away-ticker 0
  ]

end

; update movement process ranges depending on the hunting behavior of conspecifics
to set-ranges
  ; set ranges to slider values depending on whether or not the conspecific has not found food (-search) or has found food (-hunt)
  set attract-range attraction-range / 75
  set avoid-range avoidance-range / 75

  ifelse [food-found] of conspecific = false
  [
    set align-range alignment-range-search / 75
  ]
  [
    set align-range alignment-range-hunt / 75
  ]
end

; determine the vector for the attraction of a bat towards nearest conspecific
to attraction

  let pos-x-c [xcor] of conspecific                    ; pull conspecific's x & y positions
  let pos-y-c [ycor] of conspecific

  let attract-dx (pos-x-c - xcor)                      ; calculate vector x & y as difference between positions of the two bats
  let attract-dy (pos-y-c - ycor)
  set attract-vect list attract-dx attract-dy          ; x & y components added to vector list

  ifelse dist-consp > avoid-range [                    ; weighing of vector determined by relative distance away from bat and avoidance range
    set attract-strength (dist-consp - avoid-range) / (attract-range - avoid-range)  ; if further away than avoidance range adjust weight (stronger further, weaker closer)
  ]
  [
    set attract-strength 0                             ; if closer than avoidance range, set weight to min (0)
  ]

  ; modify attraction strength exponentially using the calibrated modifier
  ifelse [food-found] of conspecific = false [
    set attract-strength (1 - (1 - attract-strength)^ attract-mod-search)
  ]
  [
    set attract-strength (1 - (1 - attract-strength)^ attract-mod-hunt )
  ]
end

; determine the vector for the alignment of bat with nearest conspecific
to alignment

  ifelse distance roost 0 > 2 [                                                          ; don't trigger if too close to the roost (to avoid issues with crowding)
    let heading-consp [heading] of conspecific                                           ; pull conspecific heading...
    let speed-consp [step-length] of conspecific                                         ; ...and step length

    ; predict future position of conspecific based on these values
    let fut-x-consp ([xcor] of conspecific) + ((sin heading-consp) *  speed-consp)
    let fut-y-consp ([ycor] of conspecific) + ((cos heading-consp) *  speed-consp)

    ; pull current distance to conspecific
    let curr-dist dist-consp

    ; assess distance to future conspecific position
    let t-dx fut-x-consp - xcor
    let t-dy fut-y-consp - ycor

    let fut-dist distancexy fut-x-consp fut-y-consp

    ; check for situations when math gets tricky
    let cond-flip 0
    if fut-dist > step-length + curr-dist [ set cond-flip 1 ]                            ; no solutions, the circles are separate
    if fut-dist < abs(step-length - curr-dist)[ set cond-flip 2 ]                        ; no solutions because one circle is contained within the other
    if fut-dist = 0 and step-length = curr-dist [ set cond-flip 3 ]                      ; circles are coincident and there are an infinite number of solutions

    ; initialize alignment vector variables
    let align-dx 0
    let align-dy 0

    ; debugging
    if debug = 3 [
      let t-bat random one-of [who] of bats
      if who = t-bat [
        print (word "bat-id:          " who)
        print (word "future xy consp: " fut-x-consp " " fut-y-consp)
        print (word "curr-dist:       " curr-dist)
        print (word "t-dxy:           " t-dx " " t-dy)
        print (word "fut-dist:        " fut-dist)
        print (word "cond-flip:       " cond-flip)
        print " "
      ]
    ]

    ; if no tricky math situations
    ifelse cond-flip = 0 [
      ; estimate triangle parameters to assess which direction would result in a equal current and future distance to the conspecific
      let fut-dist-a (step-length * step-length - curr-dist * curr-dist + fut-dist * fut-dist)/(2 * fut-dist)
      let height sqrt(step-length * step-length - fut-dist-a * fut-dist-a)

      let x-mid xcor + fut-dist-a * t-dx / fut-dist
      let y-mid ycor + fut-dist-a * t-dy / fut-dist

      ; two points will be possible
      let x-new-1 x-mid + height * t-dy / fut-dist
      let x-new-2 x-mid - height * t-dy / fut-dist
      let y-new-1 y-mid - height * t-dx / fut-dist
      let y-new-2 y-mid + height * t-dx / fut-dist

      ; estimate the heading towards each
      let head-1 towardsxy x-new-1 y-new-1
      let head-2 towardsxy x-new-2 y-new-2

      ; if conspecific is hunting check then bats can maintain their current direction rather than both facing the same direction
      if [food-found] of conspecific = true [
        let diff subtract-headings heading-consp heading
        if abs diff > 90 [
          set heading-consp heading-consp - 180
        ]
      ]

      ; estimate difference in headings to two points
      let diff-head-1 abs (heading-consp - head-1)
      let diff-head-2 abs (heading-consp - head-2)


      ; use the point that is closer to the desired heading
      ifelse diff-head-1 < diff-head-2 [
        set align-dx (sin (towardsxy x-new-1 y-new-1)) * step-length
        set align-dy (cos (towardsxy x-new-1 y-new-1)) * step-length
      ]
      [
        set align-dx (sin (towardsxy x-new-2 y-new-2)) * step-length
        set align-dy (cos (towardsxy x-new-2 y-new-2)) * step-length
      ]

      ; debugging
      if debug = 3 [
        let t-bat random one-of  [who] of bats
        if who = t-bat [
          print (word "bat-id:          " who)
          print (word "CF1: fut-dist-a & height: " fut-dist-a " " height)
          print (word "CF1: xy-mid:              " x-mid " " y-mid)
          print (word "CF1: head-1&2:            " head-1 " " head-2)
          print (word "CF1: diff-head-12:        " diff-head-1 " " diff-head-2)
          print " "
        ]
      ]
    ]

    [
      ; if circles are separate just head towards the future point
      if cond-flip = 1 [
        set align-dx (sin (towardsxy fut-x-consp fut-y-consp)) * step-length
        set align-dy (cos (towardsxy fut-x-consp fut-y-consp)) * step-length
      ]
      ; when other conditions are true, stop aligning
      if cond-flip = 2 or cond-flip = 3 [
        set align-dx 0
        set align-dy 0
      ]
    ]
    set align-vect list align-dx align-dy              ; x & y components added to vector list
  ]

  [
    set align-vect list 0 0                            ; don't align if too close to roost
  ]

  ifelse dist-consp > align-range [                                                  ; weighing of vector determined by relative distance away from bat and alignment range
    set align-strength (attract-range - dist-consp) / (attract-range - align-range)  ; if further away than alignment range adjust weight (stronger closer, weaker further)
  ]
  [
    set align-strength (dist-consp - avoid-range) / (align-range - avoid-range)      ; if closer than alignment range adjust weight (weaker closer, stronger further)
  ]

  ; modify alignment strength exponentially using the calibrated modifier
  ifelse [food-found] of conspecific = false [
    set align-strength (1 - (1 - align-strength)^ align-mod-search)
  ]
  [
    set align-strength (1 - (1 - align-strength)^ align-mod-hunt)
  ]

end

; determine the vector for the avoidance of nearest conspecific if bat gets too close
to avoidance

  ifelse dist-consp < align-range [                    ; if bats are within the alignment range zone then calculate

    let pos-x-c [xcor] of conspecific                  ; pull conspecific's x & y positions
    let pos-y-c [ycor] of conspecific

    let avoid-dx (xcor - pos-x-c)                      ; calculate vector x & y as the inverse of the difference between positions of the two bats
    let avoid-dy (ycor - pos-y-c)
    set avoid-vect list avoid-dx avoid-dy              ; x & y components added to vector list

    ifelse dist-consp > avoid-range [                  ; weighing of vector determined by relative distance away from bat and avoidance range
      set avoid-strength (align-range - dist-consp) / (align-range - avoid-range) ; if further away than avoidance range adjust weight (stronger closer, weaker further)
    ]
    [ set avoid-strength 1 ]                           ; if closer than avoidance range, set to weight to max (1)
  ]
  [ set avoid-strength 0 ]                             ; if bats aren't too close, then don't use avoidance

  ; modify avoidance strength exponentially using the calibrated modifier
  ifelse [food-found] of conspecific = false [
    set avoid-strength (1 - (1 - avoid-strength)^ avoid-mod-search)
  ]
  [
    set avoid-strength (1 - (1 - avoid-strength)^ avoid-mod-hunt)
  ]
end

; determine the vector for RW and biased-RW towards foraging area
to random-walk

  let heading-rw-head heading                     ; set heading as previous heading
  let target foraging-area                        ; animals outside of the foraging region then set the foraging area center as a target
  let heading-rw-forg towards target              ; set heading towards center of foraging area

  let dist-forg 0
  let rad1 foraging-radius - 20                                                                   ; the radius where target is set to 100% current heading
  let rad2 foraging-radius + 5                                                                    ; the radius beyond where target is 100% center of the foraging region
  ifelse (((distance foraging-area - rad1) < (rad2 - rad1)) and (distance foraging-area > rad1))[ ; if not within foraging region but closer than the max distance to a side
    set dist-forg (distance foraging-area - rad1) / (rad2 - rad1)                                 ; weigh RW vector magnitude based on distance
  ]
  [
    ifelse (distance foraging-area < rad1) [ set dist-forg 0 ] [ set dist-forg 1 ]
  ]

  ; update heading based on the target
  let head heading-rw-head + ((dist-forg * 0.1) * (subtract-headings heading-rw-forg heading-rw-head))

  ; debugging
  if debug = 4 [
    let t-bat random one-of [who] of bats
    if who = t-bat [
      print (word "bat-id:     " who)
      print word "my head:     " heading-rw-head
      print word "target head: " heading-rw-forg
      print word "diff head:   " (subtract-headings heading-rw-head heading-rw-forg)
      print word "dist-forg:   " dist-forg
      print word "final head:  " head
    ]
  ]

  let ta-mod 0
  ifelse null-model? = true or conspecific = nobody  ; adjust the wiggle added to random walk when having a conspecific
  [ set ta-mod 1 ]
  [ set ta-mod 0.45 ]

  ; calculate additional wiggle to add
  set turning-angle ((0 + ( random-gamma turn-ang-alpha-search turn-ang-lambda-search ) - ( random-gamma turn-ang-alpha-search turn-ang-lambda-search )) * (180 / pi)) * ta-mod

  ; update heading and calculate RW vector
  let heading-rw head + turning-angle
  let rw-dx (sin heading-rw * step-length)            ; calculate vector if the bat were to take a step using this heading
  let rw-dy (cos heading-rw * step-length)
  set rw-vect list rw-dx rw-dy                        ; x & y components added to vector list

  ; modify RW strength using the calibrated modifier
  ifelse conspecific = nobody or [food-found] of conspecific = false
  [
    set rw-strength rw-mod-search
  ]
  [
    set rw-strength rw-mod-hunt
  ]
 end


; determine resulting movement vector and move
to bats-move

  ifelse flying-towards-food = false [

    let head-1 heading                                    ; pull current heading

    ; average weighed vectors
    let avg-x (((item 0 attract-vect) * attract-strength) + ((item 0 align-vect) * align-strength) + ((item 0 avoid-vect) * avoid-strength) + ((item 0 rw-vect) * rw-strength))
    let avg-y (((item 1 attract-vect) * attract-strength) + ((item 1 align-vect) * align-strength) + ((item 1 avoid-vect) * avoid-strength) + ((item 1 rw-vect) * rw-strength))

    facexy (xcor + avg-x) (ycor + avg-y)                  ; bats face resulting vector

    let head-2 heading                                    ; pull heading after updating

    set turn-angle-real subtract-headings head-1 head-2   ; calculate total turning angle in the step
  ]

  [
    ; if flying towards a found food cell, just face that cell
    face targ-cell
  ]

  ; move forward and check for food
  if food-found = false  [                              ; bats move forward if not on an unclaimed food cell
    repeat 3 [                                          ; split into 3 steps to avoid bats flying directly through a food cell but not stopping
      fd step-length / 3
      food-check
    ]
  ]
end

; bats which found a food cell fly around - using same process as random-walk procedure
to hunting-fly

  ifelse hunt-fly? = true
  [
    ; set one heading as current heading
    let heading-rw-head heading
    ; set other heading towards center of food cell
    let target food-cell
    let heading-rw-cell towardsxy (item 0 target + random-float 0.1 - random-float 0.1) (item 1 target + random-float 0.1 - random-float 0.1)

    let dist-cell 0
    let rad1 0.1    ; inner radius (only current heading)
    let rad2 0.5    ; outer radius (only cell center)
    ifelse (((distancexy (item 0 target)(item 1 target) - rad1) < (rad2 - rad1)) and (distancexy (item 0 target)(item 1 target) > rad1))[ ; if not within foraging region but closer than the max distance to a side
      set dist-cell (distancexy (item 0 target)(item 1 target) - rad1) / (rad2 - rad1) ; weigh RW vector magnitude based on distance
    ]
    [
      ifelse (distancexy (item 0 target)(item 1 target) < rad1) [ set dist-cell 0 ] [ set dist-cell 1 ]
    ]

    let head heading-rw-head + (dist-cell * (subtract-headings heading-rw-cell heading-rw-head))

    ; turning angles and step lengths based on hunting movements from data fit with a gamma distribution
    set turning-angle ((0 + ( random-gamma turn-ang-alpha-hunt turn-ang-lambda-hunt ) - ( random-gamma turn-ang-alpha-hunt turn-ang-lambda-hunt )) * (180 / pi))
    set heading head + turning-angle

    fd step-length  ; fly forward
  ]
  [
    ; if hunt-fly is false bats instead just wiggle around in the center of the food cell
    set heading heading + ((0 + ( random-gamma turn-ang-alpha-hunt turn-ang-lambda-hunt ) - ( random-gamma turn-ang-alpha-hunt turn-ang-lambda-hunt )) * (180 / pi))
  ]
end


; bats check if they have found food and if so they claim the food cell
to food-check

  let food-found-this-tick false

  if food-found = false [
    ; check if food is present in the cell currently in
    ifelse [food] of patch-here = true [
      ; if food is here, check if food cell has already been found by another bat
      if [found] of patch-here = false [
        set xcor pxcor
        set ycor pycor
        set food-found-this-tick true
      ]
      ; if food has already been found, check for empty cells in the same patch
      if [found] of patch-here = true and bats-per-patch? = true and null-model? = false and no-patch != no-food-cells and flying-towards-food = false [
        let pid [patch-id] of patch-here
        let nP count patches with [patch-id = pid and found = false]

        if nP > 0 [
          ; if there are empty cells, target one
          let tP one-of patches with [patch-id = pid and found = false]
          set targ-cell tP
          set flying-towards-food true
        ]
      ]
    ]

    [
      ; if food is not here
      if check-radius? = true [
        ; check in detection radius
        if any? patches in-radius prey-detection-range with [food = true] [
          let in-r-cells patches in-radius prey-detection-range with [food = true]
          ifelse any? in-r-cells with [found = false] [
            ; if there are empty cells, target one
            let tP one-of patches in-radius prey-detection-range with [food = true and found = false]
            set targ-cell tP
            set flying-towards-food true
          ]
          [
            if bats-per-patch? = true and null-model? = false and no-patch != no-food-cells and flying-towards-food = false [
              ; if there are empty cells in a detected patch, target one
              let pid [patch-id] of in-r-cells
              let pos-cell []
              ask in-r-cells [
                let myID patch-id
                if any? patches with [patch-id = myID and found = false] [
                  set pos-cell lput myID pos-cell
                ]
              ]
              if not empty? pos-cell [
                let tP one-of patches with [member? patch-id pos-cell and found = false]
                set targ-cell tP
                set flying-towards-food true
              ]
            ]
          ]
        ]
      ]
    ]
  ]

  ; update parameters when food is found
  if food-found-this-tick = true [
    set food-found true
    set conspecific nobody
    set food-cell list xcor ycor
    set color white
    ask patch-here [ set found true ]
    set flying-towards-food false
    set food-found-ticks-list lput (ticks - leave-roost-tick) food-found-ticks-list
    pu
  ]

  ; if local-enhancement? is on, bats which found food tell bats that currently have them as their conspecific where the food patch they found is
  if local-enhancement? = true and food-found-this-tick = true [ local-enhancement ]

  ; if targeting a cell, then arrive to find it occupied, reset flying-toward-food
  if patch-here = targ-cell and [found] of patch-here = true [ set flying-towards-food false ]

  ; check if 95% of bats have found food
  if precision (1 - ((count bats with [food-found = false]) / n-bats)) 2 >= 0.95 and most-done = false [
    set time-to-95-food ticks
    set most-done true
  ]
end

; local enhancement alerts bats to food if they come within 240m of a foraging bat
to local-enhancement
 ifelse food-found = false
  [
    ; bats which have encountered a hunting bat
    let me self
    ask conspecific [
      let consp-patch patch item 0 food-cell item 1 food-cell
      let pid [patch-id] of consp-patch
      ; if there are empty cells in the patch, they target one
      if any? patches with [ patch-id = pid and found = false ] [
        ask me [
          let tP one-of patches with [ patch-id = pid and found = false ]
          set targ-cell tP
          set flying-towards-food true
        ]
      ]
    ]
  ]
  [
    ; bats which have found food direct bats which see them as their conspecific to available cells if they exist
    let pid [patch-id] of patch item 0 food-cell item 1 food-cell
    ; if there are empty cells in the patch, they target one
    if any? bats with [ conspecific = myself ] and any? patches with [ patch-id = pid and found = false ] [
      ask bats with [ conspecific = myself ] [
        let tP one-of patches with [ patch-id = pid and found = false ]
        set targ-cell tP
        set flying-towards-food true
      ]
    ]
  ]
end

;;;; outputs and monitors ;;;;

to update-monitors
  set-current-plot "Network size"
  if any? bats [plot avg-netwrk]

  set-current-plot "Percent of searching bats with a conspecific"
  if any? bats with [food-found = false] [
    let perc-search-w-c count bats with [food-found = false and conspecific != nobody]/(count bats with [food-found = false])
    plot perc-search-w-c
  ]

  set-current-plot "Ticks which bats found food"
  set-plot-x-range 1 ticks
  histogram food-found-ticks-list
end

; estimate network size of bats as groups of interacting bats
to calculate-network-size
  ask bats [ set network-ID 0 ]                                                   ; reset network-ID to zero
  set networks []                                                                 ; reset networks list
  set network-unassigned-bats [who] of bats with [food-found = false]             ; create a list of searching bats to pull from

  while [length network-unassigned-bats > 0 ] [                                   ; go through the list of searching bats to ID groups
    ask one-of bats with [ who = item 0 network-unassigned-bats ] [               ; start from first position
      let ntwrk []                                                                ; initialize a group specific network list
      set network-unassigned-bats remove who network-unassigned-bats              ; remove ID from list
      if conspecific != nobody and [ food-found ] of conspecific = false [        ; check if they have a searching conspecific
        set network-ID max [network-ID] of bats + 1                               ; set their network-ID as the next available number
        let id network-ID
        if debug = 2 [print word "Network ID: " id]
        set ntwrk lput who ntwrk                                                  ; add conspecific's who to group list
        if debug = 2 [print word "Initial network: " ntwrk ]

        let cont true                                                             ; boolean to end while loop
        let consp 0                                                               ; temp variable to keep track of the who of conspecific
        set consp [ who ] of conspecific
        if debug = 2 [print word "Conspecific ID: " consp ]

        ifelse [network-ID] of bat consp != 0 [                                   ; check if conspecific is already in a network group
          if debug = 2 [ print "yes1" ]
          set network-ID ([network-ID] of conspecific)                            ; if so pull network-ID from their group
          let temp-n item (network-ID - 1) networks                               ; grab the group list
          if debug = 2 [ print temp-n ]
          set temp-n lput who temp-n                                              ; add the conspecific's who to list
          set ntwrk temp-n                                                        ; update group list
          if debug = 2 [ print ntwrk ]
          set networks replace-item (network-ID - 1) networks ntwrk               ; update full network list
          if debug = 2 [
            print network-ID
            print networks ]
        ]
        [
          while [cont = true] [                                                   ; if the conspecific is not already in a network group
            ask bat consp [
              set ntwrk lput consp ntwrk                                          ; add their who to group list
              set network-unassigned-bats remove consp network-unassigned-bats    ; remove their who from unassigned bats list
              set network-ID id                                                   ; set their network ID

              ifelse conspecific != nobody and [ food-found ] of conspecific = false [ ; if they have a searching conspecific
                if debug = 2 [ print "yes2" ]
                ifelse not member? ([ who ] of conspecific) ntwrk [               ; check if they are not already in network group
                  set consp [who] of conspecific                                  ; if not, add set the temp variable to their who and rerun loop
                  if debug = 2 [ print consp ]
                ]
                [ set cont false                                                  ; if yes, then stop
                  if debug = 2 [ print "no1" ]
                ]
              ]
              [ set cont false                                                    ; if yes, then stop
                if debug = 2 [ print "no2" ]
              ]
            ]
          ]
          if debug = 2 [
            if [ food-found ] of conspecific = false [ print "consp feeding" ]
            print ntwrk
          ]

          set networks lput ntwrk networks                                        ; update networks with the group network

          if debug = 2 [
            print network-ID
            print networks
          ]
        ]
      ]
    ]
  ]

  let single-search count bats with [ food-found = false and network-ID = 0 ]     ; add loner bats to network list as individuals
  repeat single-search [
    set networks lput [1] networks
  ]

  ; update colors
  ask bats [ifelse food-found = true [ set color white ] [ set color (min [who] of bats with [network-ID = [network-ID] of myself] + 10 )]]
end


to collect-outputs

  if ( remainder ticks  8 ) = 0 [
    ; count how many bats have a conspecific
    ifelse (count bats with [food-found = false]) > 0 [
      set perc-search-w-consp lput (precision (count bats with [food-found = false and conspecific != nobody]/(count bats with [food-found = false])) 3) perc-search-w-consp ]
    [ set perc-search-w-consp lput 1 perc-search-w-consp ]

    ; calculate percentage of cells found by bats
    set perc-cells-found lput (precision (count bats with [food-found = true] / no-food-cells) 3) perc-cells-found

    ; calculate average network size
    if length networks > 0 [
      let itm 0
      let temp-list []
      while [itm < length networks ] [
        set temp-list lput (length item itm networks) temp-list
        set itm itm + 1
      ]
      set avg-netwrk ( precision (mean temp-list) 3 )
      set average-network-size lput avg-netwrk average-network-size
    ]
  ]

  ; record network size of tracked bats
  ask bats with [ food-found = false and leave-roost-tick < ticks and member? who track-list = true ]
  [
    let tmp []
    if null-model? = false [
      ; record network size of each bat
      ifelse network-ID = 0
      [ set tmp (list ticks who 1) ]
      [ set tmp (list ticks who (count bats with [network-ID = [network-ID] of myself])) ]
      set network-sizes lput tmp network-sizes
    ]

    ; record nearest neighbor of tracked bats
    let nghbr min-one-of other turtles [distance myself]
    let dst precision ((distance nghbr) * 75) 3
    set tmp (list ticks who dst)
    set nearest-neighbor-distance lput tmp nearest-neighbor-distance
  ]


  if run-done = true [
    print "Final outputs:"
    set time-to-all-food ticks
    print word "Percent of searching bats with a conspecific per minute:    " perc-search-w-consp
    print word "Percent of cells found by bats per minute:                  " perc-cells-found
    print word "List of network sizes per bat:                              " network-sizes
    print word "Distance to nearest neighbor per bat:                       " nearest-neighbor-distance
    print word "Ticks taken for all bats to find food:                      " time-to-all-food
    print word "Ticks taken for 95% of bats to find food:                   " time-to-95-food
    set currseed floor ((behaviorspace-run-number - 1) / 2)
  ]
end

to movement-outputs
  ; collect movement tracks for calibration and evaluation

  ; for calibration:
  if dev-steps = "calibration" [
    ask bats [ if member? who track-list and food-found = false and flying-towards-food = false and leave-roost-tick < ticks [
      set diff-dist 0

      let rec true
      if prev-consp = 0 [ set rec false ]

      let convert 75
      if prev-consp != 0 [ set diff-dist ((distance prev-consp)- prev-dist-c) * convert ]

      if diff-dist != 0 [ if prev-consp-hunt = true [                          ; if calibrating for hunting conspecifics
      ; if diff-dist != 0 [ if prev-consp-hunt = false [                         ; if calibrating for searching conspecifics
        let out-list (list (precision prev-dist-c 3) (precision diff-dist 3))
        set track-outputs lput out-list track-outputs
        ]
      ]

      ifelse conspecific != nobody [
        set prev-consp conspecific
        ifelse [food-found] of conspecific = true [ set prev-consp-hunt true ] [ set prev-consp-hunt false ]
        set prev-dist-c (distance prev-consp)
      ]
      [
        set prev-consp 0
        set prev-consp-hunt false
        set prev-dist-c 0
      ]
      ]
    ]
  ]

  ; for evaluation:
  if dev-steps = "evaluation" [
    let convert 75
    ask bats [if member? who track-list and food-found = false and leave-roost-tick < ticks  [
      let out-list (list (ticks - leave-roost-tick) who (precision (xcor * convert) 3) (precision (ycor * convert) 3) )
      set track-outputs lput out-list track-outputs
      ]
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
704
505
-1
-1
6.0
1
10
1
1
1
0
0
0
1
-40
40
-40
40
0
0
1
ticks
30.0

BUTTON
14
10
78
43
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
80
10
143
43
NIL
Go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
145
10
208
43
NIL
Go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
709
261
1192
381
Network size
NIL
NIL
0.0
10.0
0.0
4.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

SLIDER
18
100
190
133
attraction-range
attraction-range
0
500
240.0
5
1
m
HORIZONTAL

SLIDER
18
209
192
242
avoidance-range
avoidance-range
0.01
200
0.01
0.01
1
m
HORIZONTAL

SLIDER
18
137
191
170
alignment-range-search
alignment-range-search
0
300
60.0
5
1
m
HORIZONTAL

SLIDER
18
62
190
95
n-bats
n-bats
80
160
160.0
1
1
NIL
HORIZONTAL

PLOT
709
141
938
261
Count bats which found food
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count bats with [food-found = true]"

SLIDER
15
259
187
292
no-patch
no-patch
1
213
213.0
1
1
NIL
HORIZONTAL

SLIDER
15
296
187
329
no-food-cells
no-food-cells
1
500
213.0
10
1
NIL
HORIZONTAL

SLIDER
18
173
191
206
alignment-range-hunt
alignment-range-hunt
0
300
120.0
5
1
m
HORIZONTAL

TEXTBOX
24
45
174
63
Bat parameters:
11
0.0
1

TEXTBOX
21
244
171
262
Map parameters:
11
0.0
1

PLOT
713
382
1194
502
Ticks which bats found food
NIL
NIL
0.0
10.0
0.0
2.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

SWITCH
14
441
186
474
null-model?
null-model?
1
1
-1000

PLOT
937
141
1192
261
Percent of searching bats with a conspecific
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ""

CHOOSER
15
376
107
421
debug
debug
0 1 2 3 4
0

TEXTBOX
17
426
167
444
For scenarios:
11
0.0
1

SWITCH
15
476
187
509
track-bats?
track-bats?
0
1
-1000

SLIDER
15
511
187
544
n-track-bats
n-track-bats
0
100
4.0
1
1
NIL
HORIZONTAL

INPUTBOX
707
10
804
70
attract-mod-search
0.006
1
0
Number

INPUTBOX
807
10
896
70
align-mod-search
0.02
1
0
Number

INPUTBOX
896
10
990
70
avoid-mod-search
0.1
1
0
Number

INPUTBOX
991
10
1075
70
rw-mod-search
0.4
1
0
Number

INPUTBOX
707
71
804
131
attract-mod-hunt
5.0E-4
1
0
Number

INPUTBOX
807
71
896
131
align-mod-hunt
3.0
1
0
Number

INPUTBOX
898
71
992
131
avoid-mod-hunt
0.0
1
0
Number

INPUTBOX
992
71
1075
131
rw-mod-hunt
0.0
1
0
Number

SWITCH
1173
10
1288
43
check-radius?
check-radius?
0
1
-1000

SWITCH
1082
10
1175
43
hunt-fly?
hunt-fly?
0
1
-1000

SWITCH
1082
42
1288
75
bats-per-patch?
bats-per-patch?
0
1
-1000

SLIDER
15
332
189
365
foraging-radius
foraging-radius
0
40
35.0
1
1
cells
HORIZONTAL

CHOOSER
107
376
199
421
dev-steps
dev-steps
"none" "calibration" "evaluation"
0

SWITCH
1082
75
1288
108
local-enhancement?
local-enhancement?
0
1
-1000

@#$#@#$#@
## ODD protocol for a bat sensory network formation model

This is an ODD protocol document (“Overview, Design concepts, and Details”), which provides a detailed description of our model presented in:
Roeleke, M., Schlägel, U. E., Gallagher, C. A., Pufelski, J., Blohm, T., Nathan, R., Toledo, S., Jeltsch, F., & C. C. Voigt. Insectivorous bats form mobile sensory networks to optimize prey localization. (Submitted.) 
See full ODD with figures, tables, & equations in the associated text.

## 1. Purpose and patterns
Searching for food patches within a sensory network may allow for animals to more efficiently locate patchily distributed and ephemeral resources, however there exists little knowledge on the potential benefits of this strategy and the conditions under which it is advantageous. The purpose of the model is to assess if and under which conditions sensory networking of insectivorous bats can lead to observed differences in the amount of time it takes bats to find food in varying environments. 
When developing the model we used patterns elucidated in the empirical portion of the manuscript (Fig. 4c & d), which document changes in distance to conspecifics associated with varying initial distances for model calibration (see ODD Element 7.8). Three additional movement patterns (time to first forage, distance from roost to forage cell, and straightness index) were used to evaluate the model (ODD Element 7.8).

## 2. Entities, state variables, and scales    
The model comprises two entities: landscape cells and individual bat agents. 

Landscape cells are square grid cells characterized by their position and whether they contain prey resources (food). If they do contain food, then cells are additionally characterized by the ID of the patch, defined as clusters of cells with food, they belong to (patch-ID) and whether they have been found by bat (found) (Table S2.1). Here food represents aggregations of insects forming static swarms which may persist in a location for some time. These could represent mating aggregations, emergence events, or groups driven by climatic conditions. The use of static food patches for the short timescales covered in the model (minutes to few hours) was supported by observations in the tracking data where individuals exhibited area-restricted search behavior in one spot for several minutes (data not shown). For simplicity, prey resources were expressed as a property of landscape cells (i.e., food present or not) rather than being explicitly described. 

Each grid cell covers 75x75 m. This cell size was selected as bats in the empirical tracking dataset tended to keep approximately this distance from each other when hunting (data not shown). The total spatial extent of the landscape is 80x80 square grid cells, covering an area of 6000x6000 m. This landscape size was selected to match the extent of the empirically studied foraging area. The model landscape has closed boundaries, i.e. was bounded at the extents and not toroidal. The food cell distribution is controlled by the parameter no-patch, which sets the total number of patches, or spatial aggregation level, of the landscape. This value can range between 1 and the number of food cells used (no-food-cells), leading to spatial aggregation levels of 100 - 0%, respectively. 

Bat agents are characterized by 17 variables related to their movement behavior and interactions with other bats (Table S2.1). Bats move through the environment using four vector-based movement behaviors: a random walk process, and three processes which are influenced by interactions with other bats (attraction, alignment, and avoidance). Each vector contains a direction and strength component. Directions point to a local x and y position (attract-vect, align-vect, avoid-vect, and rw-vect), while strengths, or weights, range between zero to one (attract-strength, align-strength, avoid-strength, and rw-strength), where one is the strongest possible value (see ODD Element 7 for details). When bats have a conspecific (conspecific), their resulting movement direction (turn-angle-real) is determined from the interacting movement processes (Fig. S2.1). Bats without conspecifics only use the random walk behavior. Their movement speed (step-length) is pulled from a distribution consistent with steps taken by bats in the empirical tracking data.

Bats search for food as they move in the landscape. When they find a food cell (food-found = “true”), they record the food’s location (food-cell) and fly in the vicinity of that cell for the remainder of the simulation period. If a bat encounters an occupied food cell in a patch which contains food cells that are unoccupied, it can target an unoccupied cell using the variables targ-cell and flying-towards-food (see food-check section for details)

Time is modelled using discrete time steps, where each step represents 8 seconds. This interval was used as it matched the sampling interval of the tracking data in the empirical portion of the study. Model runs continue until all bats have found an available food cell, generally within the span of minutes to several hours depending on the number of bats and landscape values used. 

## 3. Process overview and scheduling
Processes: To simulate the movements of interacting bats, the model proceeds using seven key processes: one related to selecting neighboring bats as conspecifics (set-conspecific), five related to calculating movement behavior (attraction, alignment, avoidance, random-walk, hunting-fly, and bats-move), and one process for identifying if a food cell has been found (food-check). The submodels for each process are described in detail in ODD Element 7. 

Bats proceed through these processes once per time step, but only run attraction, alignment, and avoidance procedures if they have not found food and have identified a conspecific in the set-conspecific step. Bats who have found food run the hunting-fly procedure once per time step.
 
Schedule: Bats in the model update their movements once per time step. Procedures all occur in the same predetermined order (Figure S3.1). All bats which have not yet found food first check for neighboring bats using the set-conspecific procedure. Then bats pull their step length from a random gamma distribution based on their foraging status (details in ODD Element 4.9). All bats which have a conspecific are first asked in a random order to calculate their attraction, alignment, avoidance, and random-walk vectors, then, in a separate step, all bats (again in a random order) are asked to move. These steps are executed separately so that bats are all on the same schedule, i.e., all calculate their movement direction, then all move (i.e., synchronous updating of their positions). Bats which have not found food and have no conspecific skip the attraction, alignment, and avoidance procedures and only run random-walk to determine their movement direction. 

Once bats have calculated their movement vectors, they then determine their resulting direction in the bats-move procedure (based on the direction and strengths of each movement vector) and move. Moving and searching for food are broken up into three steps to avoid bats passing through but missing a food cell. 

Bats which have found food do not check for conspecifics nor do they calculate any of the movement processes. They instead only fly around the cell they have found food in using the hunting-fly procedure. 

The model proceeds until all bats have found food. 

If producing outputs for calibration or evaluation, the movement-outputs procedure is executed at the beginning of the time step either once per 4 time steps (for calibration) or once every time step (for evaluation). For scenario outputs, once per every 8 time steps (~ a minute) observer outputs are calculated at the end of the time step using the collect-outputs procedure. 

## 4. Design concepts 
4.1. Basic principles
The model is built on the underlying theory that bats form sensory networks which may help them find food patches, as first proposed by Cvikel et al. (2015). The empirical portion of this study identified patterns in the movement behavior of neighboring insectivorous common noctule bats (Nyctalus noctula) which indicate the formation of sensory networks. The model was developed to quantify a potential benefit of these networks, a reduction in the length of time it takes bats to find food. 

At the agent level, the model is rooted in Boids movement dynamics (Reynolds, 1987), where bats interact using three movement processes: attraction (bats are attracted towards conspecifics), alignment (bats maintain distance from conspecifics), and avoidance (bats turn away from conspecifics). Each movement process is represented as a separate submodel (ODD Element 7). This model differs from standard Boids models in that behaviors are based only on the closest conspecific rather than all other bats within a certain radius. This was done to be able to directly use the patterns established in the empirical portion of this study (main text Fig. 4c & d) for calibrating the movement processes, where only the nearest tagged conspecific was registered. 

4.2. Emergence
The movement behavior and networking of bats emerge from bat movement decisions. The time it takes to find a food cell emerges from the landscape structure (number and placement of patches) and interactions between bats. 

4.3. Adaptation
Bats adapt their behavior through adjustments to their movement direction: bats base their direction on the distance to their conspecific (if any) and the resulting weighted headings given by each of the four movement processes. The resulting behavior is determined primarily by the movement process which is the strongest at the current distance to the conspecific. The strength of each movement process at varying distances to the conspecific was calibrated using the empirical patterns (see ODD Element 7.8 for details). 

4.4. Sensing
Bats can sense their location and the location, movement direction, and foraging behavior of other bat agents in their vicinity (240m). They can also sense whether food is present in cells within a 15m radius of their position (prey-detection-range).

4.5. Interaction
Bats directly interact with one another by influencing the movement direction of other bats. Bats only interact with one conspecific at a time in the model, forming dyads (see set-conspecific submodel below for details on how conspecifics are selected). 

In nature, bats can use echolocation calls from conspecifics to distinguish the foraging behavior of other bats and eavesdrop on local prey availability (Gager, 2019; Kalko, 1995). In the model we assume that bats can detect and distinguish the foraging behavior of conspecifics at distances of 240 m. Though this distance is 1.5 times the presumed maximum conspecific detection distance of 160 m for common noctule bats (Voigt et al., 2021), this range of 240 m was used to match the statistical analyses presented in the main text which were used for model calibration. 

Mediating interactions occur in the model through competition for food resources. Each food cell is available to only a single bat and once a food cell is occupied it becomes unavailable to all other bats for the duration of the simulation. This was done to ensure that bats could forage in food cells with sizes corresponding to the distances maintained between hunting bats in the empirical dataset. Though in nature multiple bats may be observed to forage simultaneously within the same 75x75 m area, hunting bats likely maintain distance to other bats due to factors not considered here, including acoustic jamming (Amichai et al., 2015), prey density, and competition. For this purpose we have limited one bat per food cell. 

Bats which have found food cells continue to fly around in the vicinity of the food cell (to simulate hunting activity) and continue to influence the movement behavior of other searching bats. 

4.6. Stochasticity
Bats leave the roost on a random time step between the start of the simulation and time step 37 (leave-roost-tick). This parameter staggers the emergence of bats from the roost so that all bats do not leave at the exact same time step.  

The step lengths (i.e., flying speeds) of bat agents are randomly pulled once per time step from a gamma distribution fit to the step lengths of tracked bats in the empirical dataset. Step length is based on foraging activity, such that bats pull from gamma distributions fit to bats which were exhibiting their current foraging behavior (i.e., hunting or searching). The maximum step length (3.4 cells tick-1) was set as 10% higher than the maximum empirical step length recorded (3.042 cells tick-1). Turning angles of bats in the random walk process (ODD Element 7.4.4) were also determined in this manner. 
Food patches are placed randomly within the foraging zone (in radius foraging-radius of the landscape center). Each map generation results in a unique landscape containing the specified number of patches (no-patch). 

4.7. Observation
Graphical output on the model interface shows the number of bats which have found food, the percentage of searching bats with a conspecific, the average network size, and the time it took each bat to find a food cell. 

A subset of bats can be tracked for calibrating and evaluating movement processes. For these tracked bats, the distance to their conspecific in the previous and current recorded time step (32 sec time intervals were used for calibration), difference in distances between the two recordings, and x and y positions are all observed. 
Network sizes were collected by assessing the number of bats which occur in a single network in a time step. Though bats can only have one conspecific at a time, larger networks can emerge from bats which have a different conspecific than the bat which sees them as its conspecific (e.g., bat 1 can see bat 2 as its conspecific, while bat 2 sees bat 3 as its conspecific, resulting in a network size of 3).  

For scenarios, the time it took all bats to find food, for 95% of bats to find food, and for each bat to find food (adjusted for the time step at which it left the roost) were recorded. 

## 5. Initialization
Upon initialization, a number of simulated bats are created based on the value of the input parameter n-bats. All bats are identical and are generated with all state variables initialized at 0, “false“, or “nobody”, where relevant. The one exception is the leave-roost-tick parameter which is selected randomly with a maximum number corresponding to 37 time steps ( ~5 minutes). 

Bats are generated at the location of the roost (coordinate = (0, -31.33)), which was based on the roost position in the real foraging area in the Uckermark site (Fig. S5.1). Bats do not start forming networks until after leaving the roost.

Food cells are placed randomly within the foraging area using the setup-maps submodel. Landscapes can vary in their number and location of patches (Figure S5.2). The number of food cells used (213) was selected as the mean number of 75x75 m cells occupied by bats in home range maps in the empirical dataset for days where a minimum of five bats were tracked concurrently. Landscapes are generated as follows: all patches are initialized with their food and found variables (Table S2.1) set to “false“ and their patch-id set to 0. The roost is then created and the first cell for each patch is then placed in the landscape at a minimum distance of 2 cells from other food cells. Then remaining food cells are placed one at a time in each patch, until all food cells are placed in the landscape. Each food cell must neighbor a cell in its patch but also not touch a food cell from another patch on any of its eight sides, ensuring that food patches do not touch. 

Simulation experiments are identical in setup beyond varying parameter values, e.g., running with 160 versus 80 simulated bats. 

## 6. Input data
The model does not use input data to represent time-varying processes.

## 7. Submodels

7.1. Check for conspecific (set-conspecific)
At the beginning of each time step searching bats select their conspecifics (if any) (Fig. S7.1). Bats first check if there are any other bats within the maximum range (240m). If so, all bats within this region are saved to a temporary list called conspecifics. If there are no bats within this region, then the bat sets its conspecific to “nobody” and proceeds to the random-walk submodel (Fig S7.1 left). Bats which have been networking with hunting conspecifics for too long will ignore conspecifics that are currently hunting when selecting their conspecific (see submodel check-consp-hunting below for details). To do this, bats which have identified any possible conspecifics in their area first check if they are currently moving away from hunting bats (i.e., move-away-ticker > 0). If so, then the bat looks for the closest searching bat in the conspecifics list and sets it as its conspecific (Fig S7.1 center). If there are no searching bats in the list, then the conspecific is set to “nobody”. If the conspecifics list is not empty and the bat has a move-away-ticker value of 0 (i.e., not moving away from hunting bats), then it chooses the closest bat in its radius as its conspecific (regardless of foraging behavior of the other bat)(Fig S7.1 right).  

7.2. Check if current conspecific is hunting (check-consp-hunting)
In the check-consp-hunting procedure bats keep track of how long they have had a conspecific which has been hunting. This procedure is run as a safety net for bats to disengage from hunting bats to avoid getting stuck in networks with bats which may be fixed in a location that does not contain any empty food cells. Flying around a hunting bat can help bats find empty food cells nearby when present, particularly when the patch size is large, but can also be a distraction to searching bats when no available food cells are nearby. In the model, bats will stay with a hunting conspecific for 40 time steps (~ 5 minutes), using a parameter called hunting-consp-ticker to keep track. After this timer reaches 40 time steps, bats then ignore hunting bats when selecting their conspecifics again for 40 time steps, keeping track using the parameter move-away-ticker. Move-away-ticker is used in the set-conspecific submodel (above) when selecting a conspecific from neighboring bats. 

7.3. Set ranges used (set-ranges)
This submodel is used to convert the ranges for each of the three Boids-based movement processes from distance in meters to distance in units of cell size. Additionally, as the calibration resulted in different values used for the alignment range for bats depending on whether they have a hunting or searching conspecific, this process is also used to select the correct range to use and convert that to distance in cells. 

7.4. Movement processes:
The following four submodels all relate to the movement behavior of modelled bats. Attraction, alignment, and avoidance are carried out only by bats which have a conspecific, while the random walk submodel is executed by all searching bat agents. All submodels consist of two parts: first a vector is drawn that indicates the heading identified by each process, then the strength (i.e., weight) of the vector is calculated. The directions and strengths of all vectors are taken together to determine the final heading of the bat in the time step.

7.4.1. Attraction
Bats which are attracted to other bats orient their heading in the direction of their conspecific. They do this by first pulling the x and y positions of their conspecific as temporary variables. The focal bat then draws the attraction vector in local space spanning from its position to the position of the conspecific (Fig. S7.2). 
After the vector is established, the vector strength is initially calculated based on the distance to the conspecific (with the highest value of one found at the attraction-range distance and the minimum value of zero at avoidance-range). This is calculated as:
attract–strength = dist–consp - avoidance–rangeattraction–range - avoidance–range
This initial vector strength value is then converted to an exponential scale using the calibrated modifier values, with the value used based on the current foraging activity of the conspecific, i.e., if the conspecific is not hunting then attract-mod-search is used, if it is hunting then attract-mod-hunt is used.

7.4.2. Alignment
Bats attempt to align their movement direction with that of their conspecific while trying to maintain distance with their conspecific when drawing their alignment vector. The process of calculating the alignment vector is outlined in Figure S7.4. Aligning bats first pull the heading (heading-consp) and speed (speed-consp) of their conspecific. This is used to then predict the position where the conspecific will be after it has moved (Fig. S7.3a). As bats do not adjust their headings for the current time step until the bats-move submodel, which is executed later as a separate step, this location does not necessarily represent exactly where the conspecific will be after moving, but instead depicts the point where they would be if they maintained their movement direction from the previous time step. 

Bats assess the current distance to their conspecific (curr-dist) and the distance between their current position and the predicted future position of the conspecific (fut-dist). The alignment vector can only point to one of two points, the two intersection points of a circle centered on the future position of the conspecific with a radius of curr-dist (Fig. S7.3c) and a circle centered on the focal bat with a radius of its step length (Fig. S7.3d). To locate these points, two triangles are drawn between the position of the focal bat, the future location of the conspecific, and the intersection points (Fig. S7.3e). Two lines are needed to assess the location of the intersection points, denoted as a and h in Fig. S7.3e. To estimate a (or fut-dist-a) the equation is used: 

Using this value h (or height) is calculated as:

The intersection point between these two lines (x-mid, y-mid) is then found using the current x and y position of the focal bat (xcor and ycor) as: 

Two new x and y positions are then identified as:

Finally, the heading towards each of the two points from the focal animal’s current position is calculated and the difference between these two headings is obtained. For animals which have a searching conspecific, the point is taken which results in a heading which has a smaller deviation from the conspecific’s current heading (heading-consp). If the conspecific is hunting, a check first occurs to see if the difference in the current headings between the two bats is greater than 90 degrees (meaning that they are currently heading in nearly opposing directions). If so, bats then take the point which is closer to the inverse of heading-consp. This allows bats to continue their path around a hunting bat, rather than mimicking the erratic turning angles that hunting bats exhibit. The selected new position is then used to draw the alignment vector (Fig. S7.3f).

Some conditions may exist where this math is not possible. One being that the circles are separate. This can occur when the step-length of the focal animal is too small to create a circle which intersects the circle around the future position of the conspecific. When this occurs, bats instead draw a vector straight towards the future position of the conspecific with a length of their step-length. This will not maintain distance but can minimize increases in distance. Additionally, circles can be concentric causing no possible solutions, or be coincident with an infinite number of solutions. In both of these cases the alignment vector is set to (0,0) and alignment does not occur.  

After the alignment vector is established, the vector strength is calculated based on the distance to the conspecific. Alignment differs from the other two processes in that it is strongest (a value of one) at the alignment-range and then decreases to zero in both directions, so that it has a value of zero at both the attraction- and avoidance-ranges (Fig. S7.4). The initial align-strength value is calculated as:

This weight is then converted to an exponential scale as in the Attraction submodel using the calibrated modifier values, align-mod-search for bats with searching conspecifics and align-mod-hunt for bats with hunting conspecifics, as:

7.4.3. Avoidance
Avoidance is used by bats to avoid collisions and potential echolocation jamming (Amichai et al., 2015) from neighboring bats. The avoidance vector is essentially calculated as the inverse of attraction. Bats first pull the x and y positions of their conspecific to get a vector heading towards the bat. Then this vector is flipped to head directly away from the conspecific (Fig. S7.5). 

The initial strength of the avoidance vector is then calculated using the distance from the conspecific where the highest value (1) occurs at the avoidance-range and then decreases to zero at the alignment-range. Avoidance is set to zero at ranges greater than the alignment-range. This strength is again modified exponentially, as for the attraction and alignment vectors, using the calibrated parameters avoid-mod-search (when conspecific is searching) and avoid-mod-hunt (when the conspecific is hunting), using the equation: 

An exponential scale was used for the vector strengths to calibrate the movement processes to the data in a way in which the overall magnitudes of each process was the same, but the distances at which the maximum value of one occurred changed between processes (Fig. S7.6). This approach modified only the shape of the curve rather than the range of magnitudes of the vector strengths. 

7.4.4. Random walk (random-walk)
The random walk submodel is the only movement process calculated for all searching bats in every time step. This process entirely controls the movement direction of bats when they have no conspecific or when running the null model. Both correlated and biased random walk behavior is used in the process. Bats select either the center of the foraging area or their previous heading to target when estimating their movement direction, based on their current distance to the center of the foraging area. When well within the foraging area (< foraging-radius - 20 cells) bats use their previous heading to target their random walk behavior, while when completely outside of the foraging area (> foraging-radius + five cells),  bats target a point that is 10% determined by targeting the center of the foraging area and 90% by their previous heading. For bats at a distance in the middle of these two ranges, the targeting is blended between the two directions (Fig. S7.7). The ranges for targeting behavior were selected as they allowed animals to thoroughly cover the foraging area when moving in the landscape without often running into the landscape walls. After the random walk target is determined, a random turning angle is then added to add stochasticity to the movement direction (see ODD Element 4.9). 

Bats calculating their random walk vector begin by first selecting their two targets (heading-rw-head for previous heading and heading-rw-forg for the center of the foraging area). Then the relative weight of the two targets is calculated linearly between the inner and outer range values based on the bats current distance to the center of the foraging area. The heading is then determined using the direction to the two targets and the calculated relative weight. The random contribution to the turning angle is then pulled from a gamma distribution fit to the turning angles of searching bats in the empirical dataset. To avoid the additional contribution from interactions between bats causing the realized turning angle of bats to deviate from observed turning angles, the contribution of additional stochasticity in turning angles is reduced by 55% for bats which currently have a conspecific. This value was selected as it resulted in total turning angles (determined as the change in heading between consecutive time steps) which resembled those of the empirical bats (Fig. S7.8). Bats then add the random component to the heading based on targeting and draw their random walk vector using this combined heading and their step length value. The strength of this vector is then set as a flat value using the calibrated parameters rw-mod-search (when conspecific is searching) and rw-mod-hunt (when the conspecific is hunting). 

7.5. Update heading and fly (bats-move)
In this submodel, bats calculate their resulting movement direction based on the four movement vectors and their strengths, and then move. To do this, bats calculate the changing x and y components of their resulting movement direction separately. To calculate the x component, bats multiply the x value from each of the vectors by its corresponding vector strength and then add together each of the four resulting values. The same is done for the y values and then the heading is updated using the resulting vector. The calculated heading is compared to the previous heading to determine the turning angle bats take in a time step (turn-angle-real).
 
Bats which have triggered the flying-towards-food behavior (meaning that they have encountered an occupied food cell with connected unoccupied food cells; details in the food-check procedure), skip the above calculations and instead face the location of their targeted cell. 

All bats which have not found food then fly in the direction of their updated heading at their flying speed (i.e., step-length) and check for food (food-check submodel). This step is broken into three parts to reduce the chance that a bat would fly directly through and miss a food cell, so bats take a step equal to their step-length divided by three and check for food, then repeat this process two more times. 

7.6. Hunting bat flight (hunting-fly)
Bats which have located a food cell fly around the cell to simulate hunting behavior. This is executed using the same mechanics as the random-walk submodel, only the center of the food cell is used rather than the center of the foraging area as a target and the ranges are changed. In this submodel, the target is entirely set to the center of the food cell when the bat is at a distance greater than 0.5 cells away from the center of the food cell, while when less than 0.1 cells away from the center bats use only their previous heading as their target. The gamma distributions used to pull the turning angle and step length values used here were fit to data from hunting bats in the empirical dataset. 

7.7. Check if food has been found (food-check)
In this submodel bats check their surroundings to locate a food cell. This process starts by first creating a temporary variable called food-found-this-tick and setting it to “false“. This variable will be used later to update the state variables of bats which have found food. Bats then check if the cell that they are currently in has food which has not been found by another bat (i.e., found = “false“). If so, then bats move to the center of that patch and flip the food-found-this-tick boolean to “true“. 

If there is food in the cell, but it is occupied by another bat, bats then use the patch-id to identify if there is another cell in the patch which has food and is unoccupied. This was done as bats are capable of detecting prey once they are in a food patch, so simulated bats were allowed to move within the swarms to locate a free cell within the patch, if available. If so, bats target one of these free cells (targ-cell) and set flying-towards-food boolean to “true”. 

If the cell instead does not contain food, bats check for food in a 15m radius around their position. If an unoccupied food cell is found, then bats target this food cell and set flying-towards-food boolean to “true”. If more than one unoccupied food cell is found, then the target cell is selected randomly. If an occupied food cell is found in the radius, then bats go through the above process for checking for unoccupied food cells in the patch containing the occupied food cell. 

When bats have found food they update a suite of variables: food-found is set to “true”, conspecific is set to “nobody”, they record the location of the food cell they now occupy (food-cell), set their color to white, ask the cell to set its found value to “true”, set flying-towards-food to “false”, and add the number of time steps that it took them to find food (adjusted for the time they left the roost) to the outputs list food-found-ticks-list.  

If a bat finds food, it alerts any bats which see it as their conspecific at that time step. If there are unoccupied food cells in the patch it found, then the other bat will target an unoccupied food cell in the patch and set their flying-towards-food boolean to “true”.

For bats which reach their targ-cell but find it occupied, they reset their flying-towards-food boolean to “false”.

7.8. Calibration and evaluation
Calibration: To ensure that the behavior of and interactions between modelled bats reflected the movement patterns found for the real bats in the study, we calibrated parameters controlling the vector strength for each of the four movement processes and the alignment range distance using an inverse modelling approach (Kramer-Schadt et al., 2007; Railsback & Grimm, 2019). The calibration was broken up into two steps with the movements of bats with and without hunting conspecifics being calibrated separately. Two related patterns were used for each step: 1) the relationship between initial distance and changes in distance between focal bats and their nearest conspecific (main text Fig. 4c & d), and 2) the overall shape of the density curve fit to the distance changes. A total of 10 parameters were calibrated, with the calibrated parameters for bats with searching conspecifics being: attract-mod-search, align-mod-search, avoid-mod-search, rw-mod-search, and alignment-range-search; and for bats with hunting conspecifics were: attract-mod-hunt, align-mod-hunt, avoid-mod-hunt, rw-mod-hunt, and alignment-range-hunt.
To calibrate these parameters, for each state (conspecific searching vs. hunting) 25 simulations were run using each possible parameter combination. The calibration for bats with searching conspecifics was run first and then the second calibration was run for bats with hunting conspecifics. In each simulation, the number of landscape patches (no-patches) was randomly selected, ranging from spatial aggregation levels of 0 - 100%. This was done as the total spatial aggregation of food resources in the real environment was unknown. For both calibrations 80 bats were generated and a number of these bats were tracked depending on the state, with 20 bats tracked for the first calibration and 80 tracked for the second. The number differed as a great deal more points were collected in the first calibration, due to the greater amount of time spent by bats networking with searching bats. 

Tracked bats executed the movement-outputs procedure to record their movement behavior once per every four time steps (32 seconds). These outputs were only collected when bats met certain criteria: they must not have found food yet (food-found = “false”), they must not currently be flying towards an unoccupied food cell (flying-towards-food = “false”), and they must have already left the roost (leave-roost-tick < time step). If these conditions are met, bats then calculate their change in distance (diff-dist) as the difference between their current distance to the conspecific they had in the previous recording (called prev-consp) and the distance they were to that conspecific in the previous recording (prev-dist-c). This method was taken to collect outputs in the same manner as they were calculated for bats in the empirical dataset. Bats then record prev-dist-c and diff-dist in an output list. After recording, bats check if they have a conspecific and, if so, save the identity (prev-consp), foraging behavior (prev-consp-hunt), and distance to the conspecific (prev-dist-c) at the end of the procedure. 

Each calibrated parameter was varied over four levels and all combinations of these four levels were tested. Ranges were identified through exploratory runs which narrowed the values tested (see Table S7.8.1). The calibration was run using the BehaviorSpace feature in NetLogo v6.2.0 (Shargel & Wilensky, 2002). 

The outputs were processed in R statistical software v4.0.3 (R Core Team, 2021) and pooled based on their parameter combination. The first pattern was assessed in the same manner as the empirical patterns by fitting a 3rd order polynomial model to 5000 points from each parameter combination using the lmer function in the “lme4” package (Bates et al., 2021, p. 4). Using the ggeffect function in the “ggeffects” package (Lüdecke et al., 2021), the resulting statistical model predictions were converted into a table with a row for each initial distance value between 0 and 240m at 10m increments. The deviation between the empirical and simulated statistical models was then estimated at each initial distance value and the root mean squared error (RMSE) was calculated for each parameter combination. 

For the second pattern, density values were calculated using the built-in density function in R at 100 equally spaced points between 0 and 350m for both the empirical and simulated differences in distances recorded. Again 5000 points were used per parameter combination. At each point, the deviation between the empirical and simulated density curves was determined and the overall error was again calculated as RMSE. 
The relative fit of the results from each parameter combination to each pattern was calculated as a ranking which increased with increasing RMSE, e.g., the parameter combination resulting in the tightest fit (i.e., lowest RMSE) held a rank of 1. These two ranks were added together and the parameter combination which yielded the lowest overall rank was selected. Calibration resulted in successful fits of the movement processes to the empirical patterns, with average model predictions for the first pattern falling within the 95% confidence intervals for both states at all distances. The selected parameter values can be found in Table S7.1 and resulting fits in Figure S7.9 below. 

The null model was additionally run (null-model? set to “true”) using the same simulation specifications to assess the pattern outputs for bats which are not networking (Fig. S7.10). 

Evaluation: To evaluate the calibrated model’s ability to emulate bat movement features, we compared model outputs to three population-level bat movement patterns. The three patterns we used were distributions of the: 1) euclidean distance between the bat starting point and its final point at a food cell (beeline), 2) beeline divided by the sum of euclidean distances between each 8 sec time step (straightness index), and 3) time difference between the bat leaving the roost and finding a food cell (time to first forage). It is important to note that these comparisons are independent, or secondary predictions, and did not include any additional parameterization or re-calibration. 
To compare these patterns, the x and y positions of simulated bats were collected once per time step for bats in landscapes with varying spatial aggregation levels (1, 2, 4, 8, 16, 32, 64, 128, or 213 patches). While 80 bats were simulated in the landscape (approximately the empirical colony size), only four bats were tracked in each simulation. Twenty-five simulations were run per spatial aggregation level, totalling in 100 tracked bats per level. 

The model outputs and empirical data were then analyzed following the same methodology using the “sp” package (Pebesma et al., 2021) in R statistical software v4.0.2 (R Core Team, 2021). The patterns were then visually compared to evaluate fit at each spatial aggregation level tested (Fig. S7.11).

The fit of the evaluation outputs varied depending on the tested spatial aggregation level. For the beeline pattern, most median values for the model outputs fell within the interquartile range of the empirical data. The only exception occurred when only a single food patch was generated, where the beeline value was recorded as being higher than the empirical pattern. The empirical observations showed a much wider spread of values, reaching distances substantially higher than recorded in the model (7584m), potentially due to natural variations in the location of food patches in the real landscape. 
Straightness index produced a U-shaped trend, with higher straightness values found at intermediate patch numbers and lower values exhibited at both low and high patch numbers. For this pattern, landscapes with between 2 and 64 patches resulted in median values which fell within the interquartile range of the empirical pattern. Though when assessing the overall shape of the distribution, landscapes with between 8 and 32 patches were observed to more similarly fit the shape of the empirical data, with a higher density of points occurring towards a straightness index value of 1.0. 

The model also resulted in a U-shaped trend for the time to first forage pattern, with both a higher median value and higher variability in output values found for low and high numbers of patches when compared to outputs for intermediate patch levels. The empirical observations predominantly occurred at very short time values, with the mean falling at only 8.8 minutes. The median values for landscapes containing between 4 and 128 patches all fell within the interquartile range of the empirical observations, with 32 patches most similarly fitting the overall distribution of the empirical pattern. 

Overall the model was capable of producing outputs similar to real bat observations from the tracking study, with the best fit occurring at intermediate spatial aggregation levels for all three tested patterns.

## References
Amichai, E., Blumrosen, G., & Yovel, Y. (2015). Calling louder and longer: How bats use biosonar under severe acoustic interference from other bats. Proceedings of the Royal Society B: Biological Sciences, 282(1821), 20152064. https://doi.org/10.1098/rspb.2015.2064
Bates, D., Maechler, M., Bolker, B., Walker, S., Christensen, R. H. B., Singmann, H., Dai, B., Scheipl, F., Grothendieck, G., Green, P., Fox, J., Bauer, A., & Krivitsky, P. N. (2021). lme4: Linear Mixed-Effects Models using “Eigen” and S4 (1.1-27.1) [Computer software]. https://CRAN.R-project.org/package=lme4
Boonman, A., Fenton, B., & Yovel, Y. (2019). The benefits of insect-swarm hunting to echolocating bats, and its influence on the evolution of bat echolocation signals. PLOS Computational Biology, 15(12), e1006873. https://doi.org/10.1371/journal.pcbi.1006873
Cvikel, N., Egert Berg, K., Levin, E., Hurme, E., Borissov, I., Boonman, A., Amichai, E., & Yovel, Y. (2015). Bats Aggregate to Improve Prey Search but Might Be Impaired when Their Density Becomes Too High. Current Biology, 25(2), 206–211. https://doi.org/10.1016/j.cub.2014.11.010
Gager, Y. (2019). Information transfer about food as a reason for sociality in bats. Mammal Review, 49(2), 113–120. https://doi.org/10.1111/mam.12146
Kalko, E. K. V. (1995). Insect pursuit, prey capture and echolocation in pipestirelle bats (Microchiroptera). Animal Behaviour, 50(4), 861–880. https://doi.org/10.1016/0003-3472(95)80090-5
Kramer-Schadt, S., Revilla, E., Wiegand, T., & Grimm, V. (2007). Patterns for parameters in simulation models. Ecological Modelling, 204(3), 553–556. https://doi.org/10.1016/j.ecolmodel.2007.01.018
Lüdecke, D., Aust, F., Crawley, S., & Ben-Shachar, M. S. (2021). ggeffects: Create Tidy Data Frames of Marginal Effects for “ggplot” from Model Outputs (1.1.1) [Computer software]. https://CRAN.R-project.org/package=ggeffects
Pebesma, E., Bivand, R., Rowlingson, B., Gomez-Rubio, V., Hijmans, R., Sumner, M., MacQueen, D., Lemon, J., Lindgren, F., O’Brien, J., & O’Rourke, J. (2021). sp: Classes and Methods for Spatial Data (1.4-5) [Computer software]. https://CRAN.R-project.org/package=sp
R Core Team. (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing. https://www.R-project.org/
Railsback, S. F., & Grimm, V. (2019). Agent-based and individual-based modeling: A practical introduction. Princeton university press.
Reynolds, C. W. (1987). Flocks, herds and schools: A distributed behavioral model. ACM SIGGRAPH Computer Graphics, 21(4), 25–34. https://doi.org/10.1145/37402.37406
Shargel, B., & Wilensky, U. (2002). BehaviorSpace. Northwestern University.
Voigt, C. C., Russo, D., Runkel, V., & Goerlitz, H. R. (2021). Limitations of acoustic monitoring at wind turbines to evaluate fatality risk of bats. Mammal Review. https://doi.org/10.1111/mam.12248
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="EvaluationBatMovementTracks" repetitions="25" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>not any? bats with [member? who track-list and food-found = false]</exitCondition>
    <metric>track-outputs</metric>
    <enumeratedValueSet variable="track-bats?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="forg-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-track-bats">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-f">
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-bats">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="5"/>
      <value value="10"/>
      <value value="50"/>
      <value value="100"/>
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-steps">
      <value value="&quot;evaluation&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ScenarioOutputTests" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>time-to-95-food</metric>
    <metric>food-found-ticks-list</metric>
    <enumeratedValueSet variable="track-bats?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="forg-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bats-per-patch?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consp-find-food-auto?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-f">
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-bats">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="40"/>
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="5"/>
      <value value="10"/>
      <value value="50"/>
      <value value="100"/>
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="CalibrationVectorModSearching" repetitions="25" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>track-outputs</metric>
    <enumeratedValueSet variable="rw-mod-for">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-mod-for">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-bats?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="forg-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attract-mod-for">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-track-bats">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="align-mod-for">
      <value value="5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="rw-mod-nonf" first="0.5" step="0.1" last="0.8"/>
    <enumeratedValueSet variable="alignment-range-f">
      <value value="120"/>
    </enumeratedValueSet>
    <steppedValueSet variable="attract-mod-nonf" first="0.001" step="0.005" last="0.016"/>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-nf">
      <value value="60"/>
      <value value="90"/>
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bats-per-patch?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="avoid-mod-nonf" first="0.05" step="0.05" last="0.2"/>
    <enumeratedValueSet variable="n-bats">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consp-find-food-auto?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="align-mod-nonf" first="0.01" step="0.01" last="0.04"/>
    <enumeratedValueSet variable="dev-steps">
      <value value="&quot;calibration&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="CalibrationVectorModHunting" repetitions="25" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>track-outputs</metric>
    <steppedValueSet variable="rw-mod-for" first="0.05" step="0.05" last="0.2"/>
    <steppedValueSet variable="avoid-mod-for" first="0" step="0.05" last="0.15"/>
    <enumeratedValueSet variable="track-bats?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="forg-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <steppedValueSet variable="attract-mod-for" first="0" step="5.0E-4" last="0.0015"/>
    <enumeratedValueSet variable="n-track-bats">
      <value value="4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="align-mod-for" first="1" step="2" last="7"/>
    <enumeratedValueSet variable="rw-mod-nonf">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-f">
      <value value="60"/>
      <value value="90"/>
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attract-mod-nonf">
      <value value="0.006"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-nf">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bats-per-patch?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoid-mod-nonf">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-bats">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consp-find-food-auto?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="align-mod-nonf">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dev-steps">
      <value value="&quot;calibration&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ScenarioOutputRevNULL" repetitions="500" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>time-to-95-food</metric>
    <metric>food-found-ticks-list</metric>
    <metric>network-sizes</metric>
    <metric>nearest-neighbor-distance</metric>
    <enumeratedValueSet variable="track-bats?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hunt-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bats-per-patch?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="local-enhancement?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-hunt">
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-bats">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="40"/>
      <value value="80"/>
      <value value="160"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="1"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
      <value value="16"/>
      <value value="32"/>
      <value value="64"/>
      <value value="128"/>
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ScenarioOutputRevNetworks" repetitions="500" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>time-to-95-food</metric>
    <metric>food-found-ticks-list</metric>
    <metric>network-sizes</metric>
    <metric>nearest-neighbor-distance</metric>
    <enumeratedValueSet variable="track-bats?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hunt-fly?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="check-radius?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bats-per-patch?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="local-enhancement?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="attraction-range">
      <value value="240"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alignment-range-hunt">
      <value value="120"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-food-cells">
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="debug">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-bats">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="40"/>
      <value value="80"/>
      <value value="160"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-patch">
      <value value="1"/>
      <value value="2"/>
      <value value="4"/>
      <value value="8"/>
      <value value="16"/>
      <value value="32"/>
      <value value="64"/>
      <value value="128"/>
      <value value="213"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foraging-radius">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="avoidance-range">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="null-model?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
