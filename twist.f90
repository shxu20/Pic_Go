program twist_for_lmp_structures
   implicit none
   character(len=256) :: filename_in, filename_out
   character(len=256) :: line, line_new
   integer :: ios, unit_in, unit_out
   integer :: atom_id, atom_type
   double precision, dimension(3) :: r, r_min, r_max, center
   double precision :: twist_angle, twist_angle_atom

   ! 输入输出文件名
   filename_in = 'relaxed_structure.lmp'
   filename_out = 'twisted.lmp'

   ! 提示用户输入扭转角度（角度），然后将其转换为弧度
   print *, 'Enter the twist angle in degrees:'
   read *, twist_angle
   twist_angle = twist_angle*3.141592653589793d0/180.0d0

   ! 初始化原子位置的最大最小值
   r_min = [1.0d30, 1.0d30, 1.0d30]
   r_max = [-1.0d30, -1.0d30, -1.0d30]

   ! 打开输入文件
   open (newunit=unit_in, file=filename_in, status='old')

   ! 首先，找到原子位置的最大最小值
   do
      read (unit_in, '(a)', iostat=ios) line
      if (ios < 0) exit  ! 到达文件末尾

      ! 尝试解析原子数据
      read (line, *, iostat=ios) atom_id, atom_type, r
      if (ios == 0) then
         r_min = min(r, r_min)
         r_max = max(r, r_max)
      end if
   end do

   ! 扭转中心
   center = 0.5d0*(r_min + r_max)

   ! 重置输入文件
   rewind (unit_in)

   ! 打开输出文件
   open (newunit=unit_out, file=filename_out, status='replace')

   ! 遍历所有行，对原子位置进行扭转
   do
      read (unit_in, '(a)', iostat=ios) line
      if (ios < 0) exit  ! 到达文件末尾

      ! 尝试解析原子数据
      read (line, *, iostat=ios) atom_id, atom_type, r
      if (ios == 0) then
         ! 计算旋转角度
         twist_angle_atom = twist_angle*(r(3) - center(3))/(r_max(3) - center(3))
         ! 扭转操作
         call twist(r, center, twist_angle_atom)
         ! 更新line内容
         write (line_new, '(i6,1x,i2,3(1x,f10.6))') atom_id, atom_type, r
         line = trim(line_new)
      end if

      ! 写入输出文件
      write (unit_out, '(a)') trim(line)
   end do

   ! 关闭文件
   close (unit_in)
   close (unit_out)
end program twist_for_lmp_structures

subroutine twist(r, center, angle)
   implicit none
   double precision, dimension(3), intent(inout) :: r
   double precision, dimension(3), intent(in) :: center
   double precision, intent(in) :: angle
   double precision, dimension(3, 3) :: rot
   double precision, dimension(3) :: r_rel, r_rel_rotated

   ! 计算相对于扭转中心的位置
   r_rel = r - center

   ! 计算旋转矩阵
   rot = reshape([dcos(angle), -dsin(angle), 0.d0, &
                  dsin(angle), dcos(angle), 0.d0, &
                  0.d0, 0.d0, 1.d0], [3, 3])

   ! 旋转
   r_rel_rotated = matmul(rot, r_rel)

   ! 更新位置
   r = r_rel_rotated + center
end subroutine twist
