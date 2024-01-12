module cart_methods
  implicit none
  private
  ! any exported subroutines should be included here
  ! public ...
  public reorder_matrix

contains
  ! make sure to purify all your subroutines before finalizing
  pure subroutine reorder_matrix(mat, r, c, split_col, split_val, split_point) bind(C, name="f_reorder_matrix")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: split_col
    integer(c_int), intent(in) :: r, c
    real(c_double), intent(in) :: split_val
    real(c_double), intent(inout) :: mat(r,c)
    integer(c_int), intent(out) :: split_point

    logical :: vec(r)
    integer(c_int) :: i

    vec = mat(:,split_col) > split_val

    !dir$ ivdep
    do i=1,c
      mat(:,i) = [pack(mat(:,i), .not. vec), pack(mat(:,i), vec)]
    end do
    split_point = count(vec)
  end subroutine reorder_matrix

  ! subroutine find_split(mat, r, c, num_to_check, out_column, out_value)
  !   ! I think I'm going to do this directly from C
  !   use, intrinsic :: iso_c_binding, only: c_int, c_double
  !   implicit none

  !   integer(c_int), intent(in) :: r, c, num_to_check
  !   real(c_double), intent(in) :: mat(r,c)
  !   integer(c_int), intent(out) :: out_column
  !   real(c_double), intent(out) :: out_value

  !   integer(c_int) :: to_check(c-1), i
  !   real(c_double) :: scores(num_to_check)

  !   ! determine which columns we'll check
  !   call shuffle_vec(to_check, c-1)

  !   !do i=1,num_to_check
  !   !  call find_gini_split(mat(:,to_check(i)), mat(:,c), r, scores, i, c-1)
  !   !end do
  ! end subroutine find_split


  pure subroutine gini_imp(classes, l, nclass, o_v)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: l, nclass
    integer(c_int), intent(in) :: classes(l)
    real(c_double), intent(out) :: o_v

    real(c_double) :: class_counts(nclass), sum_amt
    integer(c_int) :: i, j
    if(l == 0) then
      o_v = 1.0
      return
    end if

    ! i have a feeling this is slow
    ! do i=1, nclass
    !   class_counts(i) = 0.0+count(classes==i) ! cast to double for later
    ! end do

    ! total is always just l
    ! total = sum(class_counts)
    sum_amt = 1.0/l
    do i=1, l
      j = classes(i)
      class_counts(j) = class_counts(j) + sum_amt
    end do
    o_v = 1.0-sum(class_counts**2)
  end subroutine gini_imp

  subroutine double_gini_imp(classes, l, nclass, mask, mask_count, o_v)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: l, nclass, mask_count
    integer(c_int), intent(in) :: classes(l)
    logical, intent(in) :: mask(l)
    real(c_double), intent(out) :: o_v

    real(c_double) :: cc_left(nclass), cc_right(nclass), sum_amt
    integer(c_int) :: i, j

    cc_left(:) = 0.0
    cc_right(:) = 0.0
    sum_amt = 1.0
    do i=1, l
      j = classes(i)
      if(mask(i)) then
        cc_left(j) = cc_left(j) + 1.0
      else
        cc_right(j) = cc_right(j) + 1.0
      end if
    end do

    if(mask_count .ne. 0) cc_left = cc_left / mask_count
    if(mask_count .ne. l) cc_right = cc_right / (l-mask_count)

    ! return weighted gini impurity
    o_v = (mask_count/l)*(1.0-sum(cc_left**2)) + ((l-mask_count)/l)*(1.0-sum(cc_right**2))
  end subroutine double_gini_imp

  subroutine find_gini_split(v, response, l, nclass, o_v, o_gini_score) bind(C, name="find_gini_split_")
    ! Here I'm going to assume that scores are INTEGERS on scale 1:n
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    integer(c_int), intent(in) :: l, nclass
    integer(c_int), intent(in) :: response(l)
    real(c_double), intent(in) :: v(l)
    real(c_double), intent(out) :: o_gini_score, o_v

    integer(c_int) :: i, j, mloc
    real(c_double) :: total_gini, gains(l)
    logical :: tmpmask(l)

    ! test variables

    ! calculate the base gini impurity
    call gini_imp(response, l, nclass, total_gini)
    gains(:) = total_gini

    ! calculate gain for every possible split point
    do i=1, l
      tmpmask(:) = v <= v(i)
      j = count(tmpmask)
      if(j == l) then
        gains(i) = -1.0
      else
        call double_gini_imp(response, l, nclass, tmpmask, j, total_gini)
        gains(i) = gains(i) - total_gini
      end if
    end do
    gains = gains / l
    mloc = maxloc(gains, dim=1)
    o_v = v(mloc)
    o_gini_score = gains(mloc)

    ! trying something new
    ! rather than evaluate every split point, let's use a greedy approach
    ! I'm going to use simulated annealing to traverse the space

  end subroutine find_gini_split

end module cart_methods