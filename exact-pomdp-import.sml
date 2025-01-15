(*** General helper functions ***)

infix ++;
infix **;

fun secl x f y = f (x,y);

fun rl_sum (l1,[]) = l1
  | rl_sum ([],l2) = l2
  | rl_sum (l1,l2) = ListPair.map Real.+ (l1,l2);
fun a ++ b = rl_sum(a,b); 

fun rl_prod (l1,l2) = ListPair.map Real.* (l1,l2);

fun lst_pred f [] [] = true
  | lst_pred f [] _ = false
  | lst_pred f _ [] = false
  | lst_pred f (x::xs) (y::ys) = f(x,y) andalso lst_pred f xs ys;

fun rl_geq (l1,l2) = lst_pred Real.>= l1 l2; (* for concave pwlf *)
fun rl_leq (l1,l2) = lst_pred Real.<= l1 l2; (* for convex pwlf *)
fun rl_eq (l1,l2) = lst_pred Real.== l1 l2;
fun approx_eq (x: real, y, delta) = x < y+delta andalso x > y-delta;
fun dotprod (l1,l2) = foldl Real.+ 0.0 (rl_prod (l1,l2));
fun transp [] = []
  | transp ([]::_) = []
  | transp rows = map hd rows :: transp(map tl rows);
fun matVecProd (mat,vec) = map (fn list => dotprod(list,vec)) mat;
fun mat ** vec = matVecProd(mat,vec);

fun assoc ([], a) = raise Empty  (* for association lists *)
  | assoc ((x,y)::pairs, a) = if x=a then y else assoc(pairs, a);

fun crossSum (GamList) =
    let fun csPair (ll1,ll2) =
	    List.concat(map (fn l1 => (map (fn l2 => l1++l2) ll2)) ll1);
	fun csList (a, []) = a
	  | csList ([], x::xs) = csList(x, xs)
	  | csList (a, x::xs) = csList(csPair(a,x), xs)
    in csList([],GamList) end;

(*** Output functions ***)

fun printList l =
    print(concat((map (fn x => concat([Real.toString(x),","])) l)@["\n"]))

fun printListList (ll) =
    let 
	val stringll = map (map (fn x => concat([Real.toString(x),","]))) ll
	val add_new_line = map (fn l => concat([concat(l),"\n"])) stringll
    in print(concat(add_new_line)) end;

fun subsTilde #"~" = #"-"
  | subsTilde c = c; 

fun stringGam (Gam) =
    map (fn (label,vec) =>
	    label::(map (fn y => String.map subsTilde (Real.toString(y))) vec)) Gam;

fun writeStringLL (outfile, strll) =
    let
	val output = TextIO.openOut(outfile)
	fun writeStrL (sl) =
	    case sl of [] => TextIO.output(output, "\n")
		     | (x::xs) => (TextIO.output(output, x ^ ","); writeStrL(xs))
	fun writeStrLL (sll) =
	    case sll of
		[] => (TextIO.flushOut(output); TextIO.closeOut(output))
	      | (x::xs) => (writeStrL(x); writeStrLL(xs))
    in
	writeStrLL(strll)
    end;

(*** single POMDP step ***)

signature PROBLEM =
  sig
      val controls     : string list
      val observations : string list
      val obs_prob     : (string * (string * real list) list) list
      val cost         : (string * real list) list
      val trans_mat    : (string * real list list) list
      val term_cost    : (string * real list) list
  end;

functor PomdpStep (P: PROBLEM) =
  struct
    fun GamStep (u,z) =
	let val cardZ = Real.fromInt(length(P.observations)) 
	in map (fn (_,vec) => (map (fn x => x/cardZ) (assoc(P.cost,u))) ++ 
		(assoc(P.trans_mat,u)**rl_prod(assoc(assoc(P.obs_prob,u),z),vec)))
	end;

    fun enumObs (u,Gam) =
	map (fn vec => (u,vec)) (crossSum(map (fn z => GamStep(u,z) Gam) P.observations));
	     
    fun enumGam (Gam) = List.concat(map (fn u => enumObs(u,Gam)) P.controls);

    fun beliefStep (u,z) bstate =
	let val num = rl_prod(assoc(assoc(P.obs_prob,u),z),transp(assoc(P.trans_mat,u))**bstate)
	    val denom = (foldl op+ 0.0 num)
	in map (fn x => x/denom) num
	end;
  end;

structure Instruction : PROBLEM = 
  struct
    val controls = ["terminate", "continue"];
    val observations = ["correct", "incorrect"];
    val C = 1.0 and t = 0.3 and r = 0.6 and I = 0.1;
    val trans_mat = [("terminate",[]), ("continue",[[1.0,0.0],[t,1.0-t]])];
    val obs_prob = [("terminate",[("correct",[1.0, r]), ("incorrect",[0.0, 1.0-r])]),
		    ("continue",[("correct",[1.0, r]), ("incorrect",[0.0, 1.0-r])])];
    val cost = [("terminate",[0.0,C]), ("continue",[I,I])];
    val term_cost = [("terminate",[0.0,C])];
  end;

structure Tiger : PROBLEM =
  struct
    val controls = ["listen", "open-left", "open-right"];
    val observations = ["tiger-left", "tiger-right"];
    val trans_mat = [("listen",[[1.0,0.0],[0.0,1.0]]),("open-left",[[0.5,0.5],[0.5,0.5]]),
		     ("open-right",[[0.5,0.5],[0.5,0.5]])]; (* tiger randomly resets *)
    val obs_prob = [("listen",[("tiger-left",[0.15,0.85]),("tiger-right",[0.85,0.15])]),
		    ("open-left",[("tiger-left",[0.5,0.5]),("tiger-right",[0.5,0.5])]),
		    ("open-right",[("tiger-left",[0.5,0.5]),("tiger-right",[0.5,0.5])])];
    val cost = [("listen",[1.0,1.0]),("open-left",[100.0,~10.0]),("open-right",[~10.0,100.0])];
    val term_cost = [("",[0.0,0.0,0.0])];
  end;

structure Prune =
  struct
    datatype solution = Null | Feasible | Optimal;
    exception Solution_problem;

    (*** The dominationCheck routine for *concave* pwlf ***)
    fun dominationCheck (Gam) =
	let fun inner (gam, Gam') =
		let val (_,vec) = gam in
		    if (not (List.exists (secl true op=) (map (fn (_,x) => rl_geq(vec, x)) Gam')))
		    then
			gam::(List.filter (fn (_,x) => not (rl_geq(x, vec))) Gam')
		    else
			Gam' end;
	in foldl inner [] Gam end;

    (*** The "lexicographicMax" function for *concave* pwlf ***)
    fun lexicographicMin (gam, gam') =
	let val (_,vec) = gam and (_,vec') = gam'
	    fun findMin ([],[]) = gam  
	      | findMin ([],_) = gam
	      | findMin (_,[]) = gam'
	      | findMin (x::xs, y::ys) =
		if x<y then gam
		else if x>y then gam'
		else findMin(xs, ys)
	in
	    findMin(vec:real list, vec') end;

    (*** The bestVector routine for *concave* pwlf: new piece added! ***)
    fun bestVector (Gam, b) =
	let val delta = 0.001 (* CHECK: best vector can be sensitive to this value! *)
	    val pair_list = map (fn (label,vec) => ((label,vec), dotprod(b,vec))) Gam
	    fun pairMin (p1,p2) =
		let val (gam,v) = p1 and (gam',v') = p2
		in
		    if approx_eq(v,v',delta) then (lexicographicMin(gam, gam'),v)
		    else if Real.<(v,v') then (gam,v)
		    else (gam',v')
		end;	    
	    val (best_gam,_) = foldl pairMin (("",[]),1E10) pair_list
	in best_gam
	end;

    (*** The findRegionPoint routine: configure C file for concave/convex pwlf ***)
    fun findRegionPoint (gam, Gam) =
	let
	    open Array	 
	    (* Import C function "LP": project must now be compiled with MLton *)		 
	    val lp_solver = _import "LP_solve" public:
			    real array * int * real ref * real array * int * int ref -> int;
	    
	    (* Prepare inputs to C function call to GLPK routines *)
	    fun constraint_mat (gam, Gam) =
		let val veclen = List.length gam
		    val Gam' = List.filter (fn gam' => not(rl_eq(gam',gam))) Gam
		    val xterms = map (fn gam' => ListPair.map Real.- (gam,gam')) Gam'
		in List.map (fn x => x @ [~1.0]) xterms
		   @ [List.tabulate(veclen, fn i => 1.0)
		      @ [0.0]] end;
            (* example constraint_mat [[~0.5,0.5,~1.0],[~0.3,0.2,~1.0],[1.0,1.0,0.0]]; *)
            (* first element is padding as indexing in glpk starts at 1, not 0 *)
	    val mat_in =Array.fromList( [0.0] @ List.concat(constraint_mat(gam, Gam)) );
	    val mat_len = Array.length mat_in; (* padded with start zero for glpk *)
	    val obj_out = ref 0.0;
	    val x_out = Array.tabulate(List.length(gam)+1, fn i => 0.0);
	    val x_len =  Array.length x_out;   (* padded with start zero for glpk *)
	    val status = ref 0;
	    
	    (* Call C function *)
	    val _ = lp_solver (mat_in, mat_len, obj_out, x_out, x_len, status);

	    (* Return results *)
	    val xlist = Array.foldr op:: [] x_out;
	in case !status of 1 => raise Solution_problem  (* solution is undefined *)
			 | 2 => (Feasible, xlist)       (* solution is feasible *)
			 | 3 => (Null, xlist)           (* solution is infeasible *)
			 | 4 => (Null, xlist)           (* no feasible solution exists *)
			 | 5 => (Optimal, xlist)        (* solution is optimal *)
			 | 6 => raise Solution_problem  (* solution is unbounded *)
	end;

    (*** The Prune routine: two new pieces added! ***)
    fun Prune Gam =
	let fun inner [] Gam' = Gam'
	      | inner Gam [] =
		let val (_,vec) = hd Gam
		    val b_start = 1.0::List.tabulate(List.length(vec)-1, fn i => 0.0)
		    val gam' = bestVector(Gam, b_start)
		    val (_,vec') = gam'
		in inner (List.filter (fn (_,x) => not(rl_eq(x,vec'))) Gam)
			 ([gam']) end
	      | inner Gam Gam' = 
		let val (_,vec) = hd Gam
		    val (solution, b) = findRegionPoint(vec, map (fn (_,x) => x) Gam')
		    val gam' = bestVector(Gam, tl b)
		    val (_,vec') = gam'	 
		in if solution = Optimal andalso rl_eq(vec',(#2 (bestVector(gam'::Gam', tl b))))
		   then inner (List.filter (fn (_,x) => not(rl_eq(x,vec'))) Gam)
			      (gam'::Gam')
		   else inner (tl Gam) Gam' end;
	in inner Gam []
	end;
  end;

(*** Complete enumeration: the Monahan algorithm ***)

val Gamstart = Tiger.term_cost;

structure Eval = PomdpStep (Tiger);

fun iterate (Gam, 0) = Gam
  | iterate (Gam, n) =
    let val Gamnext = Prune.Prune(Prune.dominationCheck(Eval.enumGam Gam))
    in iterate(Gamnext, n-1) end;

val Gamout = iterate(Gamstart, 4);
writeStringLL ("output.dat", stringGam(Gamout));

