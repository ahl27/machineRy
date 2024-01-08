module cart_methods
  implicit none
  private
  ! any exported subroutines should be included here
  ! public ...
  public reorder_matrix

contains
  ! make sure to purify all your subroutines before finalizing
  subroutine reorder_matrix(mat, r, c, split_col, split_val) bind(C, name="f_reorder_matrix")
    use, intrinsic :: c_iso_binding, only: c_int, c_double
    implicit none
    integer(c_int), intent(in) :: r, c, split_col
    real(c_double), intent(in) :: split_val
    real(c_double), intent(inout) :: mat(r,c)

    logical :: vec(r)
    integer(c_int) :: i

    vec(:) = mat(,split_col) > split_val
    do i=1,c
      mat(:,i) = [pack(mat(:,i), .not. vec), pack(mat(:,i), vec)]
    end do
  end subroutine reorder_matrix
end module cart_methods